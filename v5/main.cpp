#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <csignal>
#include <condition_variable>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <new>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(ASSEMBLYCPP_USE_OPENMP)
    #include <omp.h>
#endif
#if defined(ASSEMBLYCPP_USE_MPI)
    #include <mpi.h>
#endif
#ifdef _WIN32
    #include <windows.h>
#endif

#if defined(ASSEMBLYCPP_USE_OPENMP)
    #define ASSEMBLYCPP_SEARCH_LOCAL thread_local
#else
    #define ASSEMBLYCPP_SEARCH_LOCAL
#endif

#if defined(ASSEMBLYCPP_USE_MPI)
int assemblyCppMpiRank = 0;
int assemblyCppMpiSize = 1;
#endif

#ifdef ASSEMBLYCPP_LIBRARY_BUILD
#include "assemblycpp.h"
#endif

using namespace std;
using vi = vector<int>;
using vb = vector<bool>;
using pii = pair<int, int>;
#include "activeWordMask.h"

#ifdef ASSEMBLYCPP_LIBRARY_BUILD
namespace assemblycpp::detail
{
#endif

constexpr int ceilLog2(int value)
{
    return std::bit_width(static_cast<unsigned int>(value - 1));
}

#include "globalPrimitives.h"
#ifdef ASSEMBLY_ENABLE_TELEMETRY
#include "searchTelemetry.h"
#endif
#include "ufds.h"
#include "molGraph.h"
#include "molfileParser.h"
#include "graphio.h"
#include "treeCanon.h"
#include "cyclicCanon.h"
#include "assemblyState.h"
#include "graphHashes.h"
#include "dagEnumeration.h"
#include "duplicateMatching.h"
#include "fragmentation.h"
#include "assemblyTranspositionTable.h"
#include "pathwayGenerator.h"
#include "searchContext.h"
#include "improvedBnB.h"
#include "signalHandler.h"
#include "ioflag.h"
#include "help.h"

/** Write the improved assembly indices recorded during the search. */
bool writeoutIntermediateMAs(const string &filename)
{
    ofstream ofs(filename);
    if (!ofs.is_open())
    {
        cerr << "error: could not open output file '" << filename << "'\n";
        return false;
    }

    for (size_t i = 0; i < intermediateMAs.size(); i++)
    {
        ofs << intermediateMAs[i].first << ' '
            << compensateDisjointAssemblyIndex(intermediateMAs[i].second) << '\n';
    }

    ofs.close();
    if (!ofs)
    {
        cerr << "error: could not write output file '" << filename << "'\n";
        return false;
    }
    return true;
}

bool hasMolfileExtension(const string &filename)
{
    return filename.ends_with(".mol");
}

/** Load one supported input without creating any output files. */
bool loadMoleculeInput(
    const string &input,
    molGraph &molGraphOutput,
    string &error
)
{
    const bool explicitMolfile = hasMolfileExtension(input);
    const string molfileName = explicitMolfile ? input : input + ".mol";
    ifstream molfile(molfileName);

    try
    {
        if (molfile.is_open())
        {
            if (verbose) cout << "Input: " << molfileName << '\n';
            molfileParser(molfile, molGraphOutput);
            return true;
        }
        if (!explicitMolfile)
        {
            ifstream graphFile(input);
            if (graphFile.is_open())
            {
                if (verbose) cout << "Input: " << input << '\n';
                graphio(graphFile, molGraphOutput);
                return true;
            }
            error = "input file not found: '" + input + "' (also tried '" +
                    molfileName + "')";
            return false;
        }
        error = "input file not found: '" + input + "'";
        return false;
    }
    catch (const std::exception &exception)
    {
        const string &parsedName = molfile.is_open() ? molfileName : input;
        error = "could not parse '" + parsedName + "': " + exception.what();
        return false;
    }
}

#ifndef ASSEMBLYCPP_LIBRARY_BUILD

bool isPrimaryProcess()
{
#if defined(ASSEMBLYCPP_USE_MPI)
    return assemblyCppMpiRank == 0;
#else
    return true;
#endif
}

#if defined(ASSEMBLYCPP_USE_OPENMP) || defined(ASSEMBLYCPP_USE_MPI)

struct ParallelReplicaResult
{
    int assemblyIndex = std::numeric_limits<int>::max();
    bool started = false;
    bool succeeded = false;
    bool runtimeLimitReached = false;
    bool enumerationLimitReached = false;
    string error;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    ParallelSearchWorkerTelemetry telemetry;
#endif
};

size_t parallelMinimumBonds()
{
    constexpr size_t defaultMinimum = 32;
    const char *configured = std::getenv("ASSEMBLYCPP_PARALLEL_MIN_BONDS");
    if (configured == nullptr || *configured == '\0') return defaultMinimum;

    size_t result = 0;
    const char *end = configured + std::char_traits<char>::length(configured);
    const auto parsed = std::from_chars(configured, end, result);
    if (parsed.ec != std::errc() || parsed.ptr != end) return defaultMinimum;
    return result;
}

bool configuredParallelBranchLeaseSize(size_t &leaseSize)
{
    const char *configured = std::getenv("ASSEMBLYCPP_BRANCH_LEASE_SIZE");
    // Zero is the internal sentinel for the adaptive default. A user-provided
    // value must remain strictly positive.
    uint64_t parsedLeaseSize = 0;
    bool valid = true;
    if (configured != nullptr && *configured != '\0')
    {
        const char *end = configured + std::char_traits<char>::length(configured);
        const auto parsed = std::from_chars(
            configured,
            end,
            parsedLeaseSize
        );
        valid =
            parsed.ec == std::errc() &&
            parsed.ptr == end &&
            parsedLeaseSize > 0 &&
            parsedLeaseSize <= std::numeric_limits<size_t>::max();
    }
    else if (configured != nullptr) valid = false;

#if defined(ASSEMBLYCPP_USE_MPI)
    const int localValid = valid ? 1 : 0;
    int allValid = localValid;
    MPI_Allreduce(
        &localValid,
        &allValid,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    if (allValid == 0)
    {
        if (isPrimaryProcess())
            cerr << "error: ASSEMBLYCPP_BRANCH_LEASE_SIZE must be a positive "
                    "integer or unset on every MPI rank\n";
        return false;
    }
    uint64_t minimumLeaseSize = parsedLeaseSize;
    uint64_t maximumLeaseSize = parsedLeaseSize;
    MPI_Allreduce(
        &parsedLeaseSize,
        &minimumLeaseSize,
        1,
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &parsedLeaseSize,
        &maximumLeaseSize,
        1,
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    if (minimumLeaseSize != maximumLeaseSize)
    {
        if (isPrimaryProcess())
            cerr << "error: ASSEMBLYCPP_BRANCH_LEASE_SIZE must have the same "
                    "value (or be unset) on every MPI rank\n";
        return false;
    }
    parsedLeaseSize = minimumLeaseSize;
#else
    if (!valid)
    {
        cerr << "error: ASSEMBLYCPP_BRANCH_LEASE_SIZE must be a positive "
                "integer\n";
        return false;
    }
#endif
    leaseSize = static_cast<size_t>(parsedLeaseSize);
    return true;
}

/** Choose about sixteen root leases per worker, bounded away from zero. */
size_t adaptiveParallelBranchLeaseSize(
    size_t rootJobCount,
    size_t workerCount
)
{
    constexpr size_t targetLeasesPerWorker =
        ParallelTaskScheduler::targetTasksPerWorker;
    workerCount = max<size_t>(1, workerCount);
    const size_t denominator =
        workerCount > numeric_limits<size_t>::max() / targetLeasesPerWorker
        ? numeric_limits<size_t>::max()
        : workerCount * targetLeasesPerWorker;
    if (rootJobCount == 0) return 1;
    const size_t sparseThreshold =
        workerCount > numeric_limits<size_t>::max() /
            parallelPromisingFrontierLeaseSize
        ? numeric_limits<size_t>::max()
        : workerCount * parallelPromisingFrontierLeaseSize;
    if (rootJobCount < sparseThreshold) return 1;
    const size_t guided = rootJobCount <= denominator
        ? 1
        : 1 + (rootJobCount - 1) / denominator;
    return guided < ParallelTaskScheduler::minimumTasksPerWorker
        ? parallelPromisingFrontierLeaseSize
        : guided;
}

int localParallelThreadCount()
{
#if defined(ASSEMBLYCPP_USE_OPENMP)
    omp_set_dynamic(0);
    return max(1, min(omp_get_max_threads(), omp_get_thread_limit()));
#else
    return 1;
#endif
}

bool hasMultipleParallelWorkers()
{
    const int localThreads = localParallelThreadCount();
#if defined(ASSEMBLYCPP_USE_MPI)
    return assemblyCppMpiSize > 1 || localThreads > 1;
#else
    return localThreads > 1;
#endif
}

bool parallelSearchEligible(const molGraph &graph)
{
    // Pathway and improvement-history output require one deterministic winning
    // search. Keep those modes serial until witness transfer is implemented.
    return
        !isPathway &&
        !writeIntermediateMAs &&
        runTimeMax == std::numeric_limits<unsigned long long>::max() &&
        hasMultipleParallelWorkers() &&
        static_cast<size_t>(graph.totalBonds) >= parallelMinimumBonds();
}

#if defined(ASSEMBLYCPP_USE_MPI)
bool mpiCommandLineOptionsAgree(const CommandLineArguments &arguments)
{
    constexpr size_t optionCount = 10;
    const array<unsigned long long, optionCount> local = {
        arguments.showHelp ? 1ULL : 0ULL,
        static_cast<unsigned long long>(ENUM_MAX),
        runTimeMax,
        isPathway ? 1ULL : 0ULL,
        removeHydrogens ? 1ULL : 0ULL,
        verbose ? 1ULL : 0ULL,
        disjointCompensation ? 1ULL : 0ULL,
        memTest ? 1ULL : 0ULL,
        writeIntermediateMAs ? 1ULL : 0ULL,
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        searchTelemetryEnabled ? 1ULL : 0ULL
#else
        0ULL
#endif
    };
    array<unsigned long long, optionCount> minimum = local;
    array<unsigned long long, optionCount> maximum = local;
    MPI_Allreduce(
        local.data(),
        minimum.data(),
        static_cast<int>(local.size()),
        MPI_UNSIGNED_LONG_LONG,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        local.data(),
        maximum.data(),
        static_cast<int>(local.size()),
        MPI_UNSIGNED_LONG_LONG,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    return minimum == maximum;
}

uint64_t parallelGraphFingerprint(const molGraph &graph)
{
    uint64_t result = UINT64_C(1469598103934665603);
    auto mix = [&](uint64_t value)
    {
        result ^= value;
        result *= UINT64_C(1099511628211);
    };
    mix(static_cast<uint64_t>(graph.totalBonds));
    mix(static_cast<uint64_t>(graph.mg.size()));
    for (const atom &entry : graph.mg)
    {
        mix(static_cast<uint64_t>(entry.type.size()));
        for (const unsigned char value : entry.type) mix(value);
        mix(static_cast<uint64_t>(entry.list.size()));
        for (const bond &edge : entry.list)
        {
            mix(static_cast<uint16_t>(edge.n));
            mix(static_cast<uint16_t>(edge.type));
        }
    }
    return result;
}

bool mpiGraphsAgree(const molGraph &graph)
{
    const uint64_t local = parallelGraphFingerprint(graph);
    uint64_t minimum = local;
    uint64_t maximum = local;
    MPI_Allreduce(
        &local,
        &minimum,
        1,
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &local,
        &maximum,
        1,
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    return minimum == maximum;
}

/** Rank-zero incumbent and cancellation state exposed through an RMA window. */
struct MpiDistributedSearchState
{
    int globalBest = std::numeric_limits<int>::max();
    int cancellationRequested = 0;
};

static_assert(std::is_standard_layout_v<MpiDistributedSearchState>);

/**
 * MPI_THREAD_FUNNELED-compatible root-request broker and incumbent exchange.
 *
 * Every worker consumes an atomically published process-local range. The MPI
 * main thread alone requests/refills that range, services rank-zero requests,
 * and progresses the RMA global bound and cancellation heartbeat.
 */
class MpiDistributedSearchController final:
    public distributedSearchController
{
    MPI_Win window = MPI_WIN_NULL;
    MpiDistributedSearchState exposedState;
    std::atomic<int> &processBest;
    ParallelTaskScheduler *scheduler = nullptr;
    const uint64_t rootJobCount;
    const uint64_t rootQueueSlotCount;
    const uint64_t refillJobCount;
    uint64_t nextRootQueueSlot = 0;
    std::atomic<uint64_t> localRootCursor;
    std::atomic<uint64_t> localRootEnd;
    std::atomic_bool rootsExhausted;
    std::atomic_bool progressRequested{false};
    std::atomic_bool incumbentDirty{false};
    std::chrono::steady_clock::time_point nextGlobalExchange{};
    bool searchFinished = false;

    static constexpr auto globalExchangeInterval =
        std::chrono::milliseconds(250);
    static constexpr int rootRequestTag = 19001;
    static constexpr int rootReplyTag = 19002;

    static uint64_t saturatedProduct(uint64_t left, uint64_t right) noexcept
    {
        return right != 0 && left > std::numeric_limits<uint64_t>::max() / right
            ? std::numeric_limits<uint64_t>::max()
            : left * right;
    }

    static uint64_t paddedRootSlotCount(
        uint64_t count,
        uint64_t blockWidth
    )
    {
        if (count == 0) return 0;
        blockWidth = std::max<uint64_t>(1, blockWidth);
        const uint64_t remainder = count % blockWidth;
        if (remainder == 0) return count;
        const uint64_t padding = blockWidth - remainder;
        if (count > std::numeric_limits<uint64_t>::max() - padding)
            throw std::length_error("distributed root queue exceeds capacity");
        return count + padding;
    }

    static void atomicMin(std::atomic<int> &target, int candidate) noexcept
    {
        int observed = target.load(std::memory_order_relaxed);
        while (
            candidate < observed &&
            !target.compare_exchange_weak(
                observed,
                candidate,
                std::memory_order_relaxed,
                std::memory_order_relaxed
            )
        ) {}
    }

    void notifyScheduler() noexcept;

    void exchangeGlobalState(bool includeCancellation) noexcept
    {
        const int localBest = processBest.load(std::memory_order_relaxed);
        int previousBest = localBest;
        MPI_Fetch_and_op(
            &localBest,
            &previousBest,
            MPI_INT,
            0,
            static_cast<MPI_Aint>(
                offsetof(MpiDistributedSearchState, globalBest)
            ),
            MPI_MIN,
            window
        );
        const int localCancellation = interruptionRequested() ? 1 : 0;
        int previousCancellation = localCancellation;
        if (includeCancellation)
        {
            MPI_Fetch_and_op(
                &localCancellation,
                &previousCancellation,
                MPI_INT,
                0,
                static_cast<MPI_Aint>(
                    offsetof(
                        MpiDistributedSearchState,
                        cancellationRequested
                    )
                ),
                MPI_MAX,
                window
            );
        }
        MPI_Win_flush(0, window);
        atomicMin(processBest, previousBest);
        if (
            includeCancellation &&
            (localCancellation != 0 || previousCancellation != 0)
        )
        {
            searchCancellationFlag.store(true, std::memory_order_release);
            notifyScheduler();
        }
    }

    std::array<uint64_t, 2> reserveRootRange(uint64_t requested) noexcept
    {
        const uint64_t begin = nextRootQueueSlot;
        if (interruptionRequested()) return {begin, begin};
        const uint64_t count = begin >= rootQueueSlotCount
            ? 0
            : std::min(
                std::max<uint64_t>(1, requested),
                rootQueueSlotCount - begin
            );
        nextRootQueueSlot = begin + count;
        return {begin, begin + count};
    }

    void serviceRootRequests() noexcept
    {
        if (assemblyCppMpiRank != 0) return;
        int available = 0;
        MPI_Status status;
        do
        {
            MPI_Iprobe(
                MPI_ANY_SOURCE,
                rootRequestTag,
                MPI_COMM_WORLD,
                &available,
                &status
            );
            if (available == 0) break;
            uint64_t requested = 0;
            MPI_Recv(
                &requested,
                1,
                MPI_UINT64_T,
                status.MPI_SOURCE,
                rootRequestTag,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE
            );
            const std::array<uint64_t, 2> reply = reserveRootRange(requested);
            MPI_Send(
                reply.data(),
                static_cast<int>(reply.size()),
                MPI_UINT64_T,
                status.MPI_SOURCE,
                rootReplyTag,
                MPI_COMM_WORLD
            );
        }
        while (available != 0);
    }

    void refillRootChunk() noexcept
    {
        if (
            searchFinished ||
            rootsExhausted.load(std::memory_order_acquire) ||
            interruptionRequested() ||
            localRootCursor.load(std::memory_order_acquire) <
                localRootEnd.load(std::memory_order_acquire)
        ) return;

        std::array<uint64_t, 2> range{};
        if (assemblyCppMpiRank == 0)
        {
            serviceRootRequests();
            range = reserveRootRange(refillJobCount);
        }
        else
        {
            MPI_Send(
                &refillJobCount,
                1,
                MPI_UINT64_T,
                0,
                rootRequestTag,
                MPI_COMM_WORLD
            );
            MPI_Recv(
                range.data(),
                static_cast<int>(range.size()),
                MPI_UINT64_T,
                0,
                rootReplyTag,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE
            );
        }
        if (range[0] < range[1])
        {
            localRootCursor.store(range[0], std::memory_order_relaxed);
            localRootEnd.store(range[1], std::memory_order_release);
        }
        else rootsExhausted.store(true, std::memory_order_release);
        notifyScheduler();
    }

public:
    MpiDistributedSearchController(
        size_t totalRootJobs,
        size_t leaseSize,
        int localWorkerCount,
        int globalWorkerOffset,
        int globalWorkerCount,
        std::atomic<int> &best
    ):
        processBest(best),
        rootJobCount(static_cast<uint64_t>(totalRootJobs)),
        rootQueueSlotCount(paddedRootSlotCount(
            rootJobCount,
            saturatedProduct(
                saturatedProduct(
                    static_cast<uint64_t>(leaseSize),
                    static_cast<uint64_t>(std::max(1, globalWorkerCount))
                ),
                parallelDistributedInitialLeasesPerWorker
            )
        )),
        refillJobCount(std::max<uint64_t>(
            1,
            saturatedProduct(
                static_cast<uint64_t>(leaseSize),
                static_cast<uint64_t>(std::max(1, localWorkerCount))
            )
        )),
        nextRootQueueSlot(std::min(
            rootQueueSlotCount,
            saturatedProduct(
                static_cast<uint64_t>(std::max(1, globalWorkerCount)),
                saturatedProduct(
                    static_cast<uint64_t>(leaseSize),
                    parallelDistributedInitialLeasesPerWorker
                )
            )
        )),
        localRootCursor(std::min(
            rootQueueSlotCount,
            saturatedProduct(
                static_cast<uint64_t>(std::max(0, globalWorkerOffset)),
                saturatedProduct(
                    static_cast<uint64_t>(leaseSize),
                    parallelDistributedInitialLeasesPerWorker
                )
            )
        )),
        localRootEnd(std::min(
            rootQueueSlotCount,
            saturatedProduct(
                static_cast<uint64_t>(
                    std::max(0, globalWorkerOffset) +
                    std::max(1, localWorkerCount)
                ),
                saturatedProduct(
                    static_cast<uint64_t>(leaseSize),
                    parallelDistributedInitialLeasesPerWorker
                )
            )
        )),
        rootsExhausted(
            saturatedProduct(
                static_cast<uint64_t>(std::max(1, globalWorkerCount)),
                saturatedProduct(
                    static_cast<uint64_t>(leaseSize),
                    parallelDistributedInitialLeasesPerWorker
                )
            ) >= rootQueueSlotCount
        )
    {
        if (totalRootJobs > std::numeric_limits<uint64_t>::max())
            throw std::length_error("distributed root queue exceeds capacity");

        if (assemblyCppMpiRank == 0)
        {
            exposedState.globalBest = best.load(std::memory_order_relaxed);
        }
        MPI_Win_create(
            assemblyCppMpiRank == 0 ? &exposedState : nullptr,
            assemblyCppMpiRank == 0
                ? static_cast<MPI_Aint>(sizeof(exposedState))
                : 0,
            1,
            MPI_INFO_NULL,
            MPI_COMM_WORLD,
            &window
        );
        MPI_Win_lock_all(0, window);
        if (assemblyCppMpiRank == 0) MPI_Win_sync(window);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    ~MpiDistributedSearchController() override
    {
        if (window == MPI_WIN_NULL) return;
        MPI_Win_flush_all(window);
        MPI_Win_unlock_all(window);
        MPI_Win_free(&window);
    }

    MpiDistributedSearchController(
        const MpiDistributedSearchController &
    ) = delete;
    MpiDistributedSearchController &operator=(
        const MpiDistributedSearchController &
    ) = delete;

    void attachScheduler(ParallelTaskScheduler &value) noexcept
    {
        scheduler = &value;
    }

    distributedRootAvailability claimRootLease(
        size_t requestedLeaseSize,
        size_t &begin,
        size_t &end
    ) noexcept override
    {
        uint64_t cursor = localRootCursor.load(std::memory_order_relaxed);
        while (true)
        {
            const uint64_t limit = localRootEnd.load(std::memory_order_acquire);
            if (cursor >= limit)
            {
                progressRequested.store(true, std::memory_order_release);
                if (rootsExhausted.load(std::memory_order_acquire))
                {
                    // The progress owner publishes cursor/end before marking
                    // the final chunk exhausted. Re-read both after acquiring
                    // that marker so a stale pre-publication limit cannot make
                    // scheduler completion sticky and lose the tail chunk.
                    cursor = localRootCursor.load(std::memory_order_acquire);
                    if (
                        cursor <
                        localRootEnd.load(std::memory_order_acquire)
                    ) continue;
                    return distributedRootAvailability::complete;
                }
                return distributedRootAvailability::wait;
            }
            const uint64_t requested = std::max<uint64_t>(
                1,
                static_cast<uint64_t>(requestedLeaseSize)
            );
            const uint64_t claimedEnd = cursor + std::min(
                requested,
                limit - cursor
            );
            if (localRootCursor.compare_exchange_weak(
                cursor,
                claimedEnd,
                std::memory_order_relaxed,
                std::memory_order_relaxed
            ))
            {
                begin = static_cast<size_t>(cursor);
                end = static_cast<size_t>(claimedEnd);
                return distributedRootAvailability::lease;
            }
        }
    }

    void requestProgress() noexcept override
    {
        progressRequested.store(true, std::memory_order_release);
    }

    void publishIncumbent(int) noexcept override
    {
        incumbentDirty.store(true, std::memory_order_release);
    }

    void progress() noexcept override
    {
        const bool requested = progressRequested.exchange(
            false,
            std::memory_order_acq_rel
        );
        const bool dirty = incumbentDirty.exchange(
            false,
            std::memory_order_acq_rel
        );
        const auto now = std::chrono::steady_clock::now();
        const bool heartbeatDue = now >= nextGlobalExchange;
        if (dirty || heartbeatDue)
        {
            exchangeGlobalState(heartbeatDue);
            if (heartbeatDue)
                nextGlobalExchange = now + globalExchangeInterval;
        }
        serviceRootRequests();
        if (
            requested ||
            localRootCursor.load(std::memory_order_acquire) >=
                localRootEnd.load(std::memory_order_acquire)
        ) refillRootChunk();
    }

    /** Keep the RMA window alive until every rank has finished local work. */
    void waitForGlobalCompletion() noexcept
    {
        searchFinished = true;
        exchangeGlobalState(true);
        // Local quiescence is permanent because descendant tasks never cross
        // rank boundaries. Keep an early rank outside a blocking collective so
        // its MPI progress calls can still service passive-target operations
        // from ranks that are finishing long roots.
        MPI_Request rendezvous = MPI_REQUEST_NULL;
        MPI_Ibarrier(MPI_COMM_WORLD, &rendezvous);
        int complete = 0;
        while (complete == 0)
        {
            progress();
            MPI_Test(&rendezvous, &complete, MPI_STATUS_IGNORE);
            if (complete == 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }
};

void MpiDistributedSearchController::notifyScheduler() noexcept
{
    if (scheduler == nullptr) return;
    if (interruptionRequested()) scheduler->cancel();
    else scheduler->notifyRootProgress();
}
#endif

bool runParallelSearch(molGraph &graph, ofstream &output)
{
    size_t branchLeaseSize = 1;
    if (!configuredParallelBranchLeaseSize(branchLeaseSize)) return false;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    const uint64_t parallelStartedNanoseconds = searchTelemetryEnabled
        ? searchTelemetryWallNanoseconds()
        : 0;
#endif
    const int localThreads = localParallelThreadCount();
    int globalWorkerCount = localThreads;
    int globalWorkerOffset = 0;
#if defined(ASSEMBLYCPP_USE_MPI)
    MPI_Allreduce(
        &localThreads,
        &globalWorkerCount,
        1,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    MPI_Exscan(
        &localThreads,
        &globalWorkerOffset,
        1,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    if (assemblyCppMpiRank == 0) globalWorkerOffset = 0;
#endif

    SearchContext searchContext;
    string preparationError;
    int preparationSucceeded = 1;
    try
    {
        prepareParallelSearchContext(
            graph,
            searchContext,
            static_cast<size_t>(localThreads)
        );
    }
    catch (const std::exception &exception)
    {
        preparationError = exception.what();
        preparationSucceeded = 0;
    }
    catch (...)
    {
        preparationError = "unknown parallel root preparation failure";
        preparationSucceeded = 0;
    }
#if defined(ASSEMBLYCPP_USE_MPI)
    int allPreparationsSucceeded = preparationSucceeded;
    MPI_Allreduce(
        &preparationSucceeded,
        &allPreparationsSucceeded,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    preparationSucceeded = allPreparationsSucceeded;
#endif
    if (preparationSucceeded == 0)
    {
        if (isPrimaryProcess())
        {
            if (preparationError.empty())
                preparationError = "another MPI rank could not prepare the root";
            cerr << "error: parallel root preparation failed: "
                 << preparationError << '\n';
        }
        return false;
    }

#if defined(ASSEMBLYCPP_USE_MPI)
    const bool distributedRootQueue = assemblyCppMpiSize > 1;
#else
    constexpr bool distributedRootQueue = false;
#endif
    const bool adaptiveBranchLeases = branchLeaseSize == 0;
    if (adaptiveBranchLeases)
    {
        branchLeaseSize = distributedRootQueue
            ? 1
            : adaptiveParallelBranchLeaseSize(
                searchContext.rootJobs.size(),
                static_cast<size_t>(max(1, globalWorkerCount))
            );
    }
    if (distributedRootQueue)
    {
        // A lease larger than the complete frontier cannot expose additional
        // work. Capping it also keeps the padded, striped queue arithmetic
        // bounded for deliberately extreme experimental overrides.
        branchLeaseSize = std::min(
            branchLeaseSize,
            std::max<size_t>(1, searchContext.rootJobs.size())
        );
    }

    std::atomic<int> processBest(searchContext.rootAssemblyIndex);
    std::atomic_bool warmStartReady(false);
#if defined(ASSEMBLYCPP_USE_MPI)
    std::optional<MpiDistributedSearchController> distributedSearchStorage;
    if (distributedRootQueue)
    {
        distributedSearchStorage.emplace(
            searchContext.rootJobs.size(),
            branchLeaseSize,
            localThreads,
            globalWorkerOffset,
            globalWorkerCount,
            processBest
        );
    }
    const size_t rankPartitionIndex = static_cast<size_t>(assemblyCppMpiRank);
    const size_t rankPartitionCount = static_cast<size_t>(assemblyCppMpiSize);
#else
    constexpr size_t rankPartitionIndex = 0;
    constexpr size_t rankPartitionCount = 1;
#endif
    unique_ptr<ParallelTaskScheduler> taskSchedulerStorage;
    vector<ParallelReplicaResult> replicas;
    int schedulerPreparationSucceeded = 1;
    try
    {
        taskSchedulerStorage = make_unique<ParallelTaskScheduler>(
            searchContext.rootJobs.size(),
            rankPartitionIndex,
            rankPartitionCount,
            static_cast<size_t>(localThreads),
            branchLeaseSize,
            adaptiveBranchLeases,
            distributedRootQueue,
            static_cast<size_t>(max(1, globalWorkerCount)) *
                parallelDistributedInitialLeasesPerWorker
        );
        replicas.resize(static_cast<size_t>(localThreads));
    }
    catch (const std::exception &exception)
    {
        preparationError = exception.what();
        schedulerPreparationSucceeded = 0;
    }
    catch (...)
    {
        preparationError = "unknown parallel scheduler preparation failure";
        schedulerPreparationSucceeded = 0;
    }
#if defined(ASSEMBLYCPP_USE_MPI)
    int allSchedulersPrepared = schedulerPreparationSucceeded;
    MPI_Allreduce(
        &schedulerPreparationSucceeded,
        &allSchedulersPrepared,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    schedulerPreparationSucceeded = allSchedulersPrepared;
#endif
    if (schedulerPreparationSucceeded == 0)
    {
        if (isPrimaryProcess())
        {
            if (preparationError.empty())
                preparationError =
                    "another MPI rank could not prepare the scheduler";
            cerr << "error: parallel scheduler preparation failed: "
                 << preparationError << '\n';
        }
        return false;
    }
    ParallelTaskScheduler &taskScheduler = *taskSchedulerStorage;
#if defined(ASSEMBLYCPP_USE_MPI)
    if (distributedSearchStorage.has_value())
        distributedSearchStorage->attachScheduler(taskScheduler);
#endif

    auto runReplica = [&](int threadIndex)
    {
        ParallelReplicaResult &result = replicas[threadIndex];
        result.started = true;
#if defined(ASSEMBLYCPP_USE_MPI)
        searchRankPartitionIndex = static_cast<size_t>(assemblyCppMpiRank);
        searchRankPartitionCount = static_cast<size_t>(assemblyCppMpiSize);
#else
        searchRankPartitionIndex = 0;
        searchRankPartitionCount = 1;
#endif
        searchBranchLeaseSize = branchLeaseSize;
#if defined(ASSEMBLYCPP_USE_MPI)
        activeDistributedSearch = distributedSearchStorage.has_value()
            ? std::addressof(*distributedSearchStorage)
            : nullptr;
        ownsDistributedSearchProgress =
            activeDistributedSearch != nullptr && threadIndex == 0;
#else
        activeDistributedSearch = nullptr;
        ownsDistributedSearchProgress = false;
#endif
        distributedSearchProgressCountdown = 0;
        parallelTaskScheduler = taskScheduler.transferableTasksEnabled()
            ? &taskScheduler
            : nullptr;
        sharedAssemblyIndex = &processBest;
        suppressSearchOutput = true;
        searchDepthTwoTasksSpawned = 0;
        searchDepthTwoTasksExecuted = 0;
        searchDeeperTasksSpawned = 0;
        searchDeeperTasksExecuted = 0;
        searchTaskStealAttempts = 0;
        searchTaskSteals = 0;
        searchLocalTaskExecutions = 0;
        searchSchedulerIdleWaits = 0;
        searchSchedulerIdleNanoseconds = 0;
        searchDeepRefillActivations = 0;
        searchTaskQueueHighWatermark = 0;
        searchMaximumTaskDepthExecuted = 0;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        searchProactiveTailRefills = 0;
#endif
        searchWarmStartBranches = 0;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        uint64_t workerStartedNanoseconds = 0;
        if (searchTelemetryEnabled)
        {
            if (threadIndex == 0)
            {
                // Worker zero owns the setup telemetry already collected by
                // this process's one root enumeration.
                workerStartedNanoseconds = parallelStartedNanoseconds;
                searchTelemetry.collectPhaseMemory = false;
            }
            else
            {
                workerStartedNanoseconds = searchTelemetryWallNanoseconds();
                // Process-wide /proc peak resets cannot safely describe
                // concurrent workers. Siblings collect counters and time only.
                resetSearchTelemetry(false);
            }
        }
#endif
        try
        {
            configureParallelWorker(
                searchContext,
                static_cast<size_t>(threadIndex)
            );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            if (searchTelemetryEnabled && threadIndex != 0)
                setSearchTelemetryPhase(SearchTelemetryPhase::assemblySearch);
#endif
            {
                // Construct and destroy every EdgeMask-owning object on this
                // worker; only the primitive job index crosses the queue.
                WorkerContext worker(
                    searchContext,
                    static_cast<size_t>(threadIndex)
                );
                if (threadIndex == 0)
                {
                    warmStartParallelIncumbent(
                        searchContext,
                        taskScheduler.warmStartRootJob(),
                        worker
                    );
                    if (activeDistributedSearch != nullptr)
                        activeDistributedSearch->progress();
                    warmStartReady.store(true, std::memory_order_release);
                    warmStartReady.notify_all();
                }
                else warmStartReady.wait(false, std::memory_order_acquire);
                if (searchContext.sharedStates != nullptr)
                {
                    runParallelRootJobs<true>(
                        searchContext,
                        worker,
                        taskScheduler
                    );
                }
                else
                {
                    runParallelRootJobs<false>(
                        searchContext,
                        worker,
                        taskScheduler
                    );
                }
                result.assemblyIndex = compensateDisjointAssemblyIndex(
                    worker.assemblyIndex
                );
            }
            clearParallelWorkerMasks();
            result.succeeded = true;
            result.runtimeLimitReached = runtimeLimitReached;
            result.enumerationLimitReached = enumerationLimitReached;
        }
        catch (const std::exception &exception)
        {
            searchCancellationFlag.store(true);
            taskScheduler.cancel();
            if (threadIndex == 0)
            {
                warmStartReady.store(true, std::memory_order_release);
                warmStartReady.notify_all();
            }
            clearParallelWorkerMasks();
            result.error = exception.what();
        }
        catch (...)
        {
            searchCancellationFlag.store(true);
            taskScheduler.cancel();
            if (threadIndex == 0)
            {
                warmStartReady.store(true, std::memory_order_release);
                warmStartReady.notify_all();
            }
            clearParallelWorkerMasks();
            result.error = "unknown parallel worker failure";
        }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled)
        {
            finaliseSearchTelemetry();
            const uint64_t workerElapsedNanoseconds =
                telemetryNanosecondDifference(
                    workerStartedNanoseconds,
                    searchTelemetryWallNanoseconds()
                );
            result.telemetry = captureParallelSearchWorkerTelemetry(
#if defined(ASSEMBLYCPP_USE_MPI)
                static_cast<uint64_t>(assemblyCppMpiRank),
#else
                0,
#endif
                static_cast<uint64_t>(threadIndex),
                static_cast<uint64_t>(globalWorkerOffset + threadIndex),
                static_cast<uint64_t>(searchRankPartitionIndex),
                static_cast<uint64_t>(searchRankPartitionCount),
                static_cast<uint64_t>(searchRootBranchOrdinal),
                static_cast<uint64_t>(searchBranchLeaseCount),
                static_cast<uint64_t>(searchBranchAssignmentCount),
                static_cast<uint64_t>(searchDepthTwoTasksSpawned),
                static_cast<uint64_t>(searchDepthTwoTasksExecuted),
                static_cast<uint64_t>(searchDeeperTasksSpawned),
                static_cast<uint64_t>(searchDeeperTasksExecuted),
                static_cast<uint64_t>(searchTaskStealAttempts),
                static_cast<uint64_t>(searchTaskSteals),
                static_cast<uint64_t>(searchLocalTaskExecutions),
                static_cast<uint64_t>(searchSchedulerIdleWaits),
                static_cast<uint64_t>(searchSchedulerIdleNanoseconds),
                static_cast<uint64_t>(searchDeepRefillActivations),
                static_cast<uint64_t>(searchTaskQueueHighWatermark),
                static_cast<uint64_t>(searchMaximumTaskDepthExecuted),
                static_cast<uint64_t>(searchProactiveTailRefills),
                static_cast<uint64_t>(searchWarmStartBranches),
                workerElapsedNanoseconds
            );
        }
#endif
        parallelTaskScheduler = nullptr;
        sharedAssemblyIndex = nullptr;
        activeDistributedSearch = nullptr;
        ownsDistributedSearchProgress = false;
        distributedSearchProgressCountdown = 0;
        suppressSearchOutput = false;
        searchRankPartitionIndex = 0;
        searchRankPartitionCount = 1;
        searchRootBranchOrdinal = 0;
        searchBranchLeaseSize = 1;
        searchBranchLeaseCount = 0;
        searchBranchAssignmentCount = 0;
        searchDepthTwoTasksSpawned = 0;
        searchDepthTwoTasksExecuted = 0;
        searchDeeperTasksSpawned = 0;
        searchDeeperTasksExecuted = 0;
        searchTaskStealAttempts = 0;
        searchTaskSteals = 0;
        searchLocalTaskExecutions = 0;
        searchSchedulerIdleWaits = 0;
        searchSchedulerIdleNanoseconds = 0;
        searchDeepRefillActivations = 0;
        searchTaskQueueHighWatermark = 0;
        searchMaximumTaskDepthExecuted = 0;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        searchProactiveTailRefills = 0;
#endif
        searchWarmStartBranches = 0;
    };

#if defined(ASSEMBLYCPP_USE_OPENMP)
    #pragma omp parallel num_threads(localThreads)
    {
        runReplica(omp_get_thread_num());
    }
#else
    runReplica(0);
#endif

#if defined(ASSEMBLYCPP_USE_MPI)
    if (distributedSearchStorage.has_value())
        distributedSearchStorage->waitForGlobalCompletion();
#endif

#ifdef ASSEMBLY_ENABLE_TELEMETRY
    uint64_t globalParallelElapsedNanoseconds = 0;
    if (searchTelemetryEnabled)
    {
        const uint64_t localParallelElapsedNanoseconds =
            telemetryNanosecondDifference(
                parallelStartedNanoseconds,
                searchTelemetryWallNanoseconds()
            );
        globalParallelElapsedNanoseconds = localParallelElapsedNanoseconds;
#if defined(ASSEMBLYCPP_USE_MPI)
        MPI_Allreduce(
            &localParallelElapsedNanoseconds,
            &globalParallelElapsedNanoseconds,
            1,
            MPI_UINT64_T,
            MPI_MAX,
            MPI_COMM_WORLD
        );
#endif
    }
#endif

    int localAssemblyIndex = compensateDisjointAssemblyIndex(
        processBest.load(std::memory_order_relaxed)
    );
    int localSucceeded = 1;
    int localRuntimeLimit = 0;
    int localEnumerationLimit = 0;
    for (const ParallelReplicaResult &replica : replicas)
    {
        localAssemblyIndex = min(localAssemblyIndex, replica.assemblyIndex);
        localSucceeded &= replica.started && replica.succeeded ? 1 : 0;
        localRuntimeLimit |= replica.runtimeLimitReached ? 1 : 0;
        localEnumerationLimit |= replica.enumerationLimitReached ? 1 : 0;
    }

    int globalAssemblyIndex = localAssemblyIndex;
    int globalSucceeded = localSucceeded;
    int globalRuntimeLimit = localRuntimeLimit;
    int globalEnumerationLimit = localEnumerationLimit;
    int globalUserInterrupt = receivedUserInterrupt() ? 1 : 0;
#if defined(ASSEMBLYCPP_USE_MPI)
    MPI_Allreduce(
        &localAssemblyIndex,
        &globalAssemblyIndex,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &localSucceeded,
        &globalSucceeded,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &localRuntimeLimit,
        &globalRuntimeLimit,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &localEnumerationLimit,
        &globalEnumerationLimit,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    const int localUserInterrupt = globalUserInterrupt;
    MPI_Allreduce(
        &localUserInterrupt,
        &globalUserInterrupt,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
#endif

#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled)
    {
        if (
            searchContext.sharedStates != nullptr &&
            !replicas.empty()
        )
        {
            const sharedAssemblyTranspositionTable::statistics stats =
                searchContext.sharedStates->stats();
            SharedAssemblyCacheTelemetry &telemetry =
                replicas.front().telemetry.sharedAssemblyCache;
            telemetry.tableCount = 1;
            telemetry.hits = stats.hitCount;
            telemetry.misses = stats.missCount;
            telemetry.collisionChainSteps = stats.collisionChainSteps;
            telemetry.allocatedBytes = stats.allocatedBytes;
            telemetry.lockAcquisitions = stats.lockAcquisitionCount;
            telemetry.lockWaits = stats.lockWaitCount;
            telemetry.lockWaitNanoseconds = stats.lockWaitNanoseconds;
        }

        vector<ParallelSearchWorkerTelemetry> localWorkerTelemetry;
        localWorkerTelemetry.reserve(replicas.size());
        for (const ParallelReplicaResult &replica : replicas)
            localWorkerTelemetry.push_back(replica.telemetry);

        vector<ParallelSearchWorkerTelemetry> gatheredWorkerTelemetry;
#if defined(ASSEMBLYCPP_USE_MPI)
        const int localTelemetryBytes = static_cast<int>(
            localWorkerTelemetry.size() *
            sizeof(ParallelSearchWorkerTelemetry)
        );
        vector<int> telemetryBytesPerRank;
        vector<int> telemetryDisplacements;
        if (isPrimaryProcess())
        {
            telemetryBytesPerRank.resize(
                static_cast<size_t>(assemblyCppMpiSize)
            );
            telemetryDisplacements.resize(
                static_cast<size_t>(assemblyCppMpiSize)
            );
        }
        MPI_Gather(
            &localTelemetryBytes,
            1,
            MPI_INT,
            isPrimaryProcess() ? telemetryBytesPerRank.data() : nullptr,
            1,
            MPI_INT,
            0,
            MPI_COMM_WORLD
        );
        if (isPrimaryProcess())
        {
            int gatheredBytes = 0;
            for (int rank = 0; rank < assemblyCppMpiSize; ++rank)
            {
                telemetryDisplacements[rank] = gatheredBytes;
                gatheredBytes += telemetryBytesPerRank[rank];
            }
            gatheredWorkerTelemetry.resize(
                static_cast<size_t>(gatheredBytes) /
                sizeof(ParallelSearchWorkerTelemetry)
            );
        }
        MPI_Gatherv(
            localWorkerTelemetry.data(),
            localTelemetryBytes,
            MPI_BYTE,
            isPrimaryProcess() ? gatheredWorkerTelemetry.data() : nullptr,
            isPrimaryProcess() ? telemetryBytesPerRank.data() : nullptr,
            isPrimaryProcess() ? telemetryDisplacements.data() : nullptr,
            MPI_BYTE,
            0,
            MPI_COMM_WORLD
        );
#else
        gatheredWorkerTelemetry = std::move(localWorkerTelemetry);
#endif
        if (isPrimaryProcess())
        {
#if defined(ASSEMBLYCPP_USE_MPI) && defined(ASSEMBLYCPP_USE_OPENMP)
            constexpr const char *parallelMode = "hybrid";
#elif defined(ASSEMBLYCPP_USE_MPI)
            constexpr const char *parallelMode = "mpi";
#else
            constexpr const char *parallelMode = "openmp";
#endif
#if defined(ASSEMBLYCPP_USE_MPI)
            constexpr const char *aggregationScope = "all_mpi_ranks";
            const uint64_t rankCount =
                static_cast<uint64_t>(assemblyCppMpiSize);
#else
            constexpr const char *aggregationScope = "process";
            constexpr uint64_t rankCount = 1;
#endif
            configureParallelSearchTelemetry(
                parallelMode,
                aggregationScope,
                rankCount,
                static_cast<uint64_t>(branchLeaseSize),
                globalParallelElapsedNanoseconds,
                globalSucceeded != 0 &&
                    globalRuntimeLimit == 0 &&
                    globalEnumerationLimit == 0 &&
                    globalUserInterrupt == 0,
                std::move(gatheredWorkerTelemetry)
            );
        }
    }
#endif

    if (globalUserInterrupt != 0)
    {
#ifdef _WIN32
        userInterruptReceived.store(true);
#else
        userInterruptReceived = 1;
#endif
    }

    lastCalculatedAssemblyIndex = globalAssemblyIndex;
    runtimeLimitReached = globalRuntimeLimit != 0;
    enumerationLimitReached = globalEnumerationLimit != 0;
    if (isPrimaryProcess() && globalSucceeded != 0)
    {
        output << globalAssemblyIndex << '\n';
        if (runtimeLimitReached) output << "status: runtime limit reached\n";
        if (enumerationLimitReached)
            output << "status: enumeration limit reached\n";
    }
    if (isPrimaryProcess())
    {
        for (const ParallelReplicaResult &replica : replicas)
        {
            if (!replica.error.empty())
                cerr << "error: parallel worker failed: " << replica.error << '\n';
        }
        if (globalSucceeded == 0)
            cerr << "error: one or more parallel workers did not complete\n";
    }
    return globalSucceeded != 0;
}

#endif

bool runConfiguredSearch(molGraph &graph, ofstream &output)
{
#if defined(ASSEMBLYCPP_USE_MPI)
    if (!mpiGraphsAgree(graph))
    {
        if (isPrimaryProcess())
            cerr << "error: MPI ranks loaded different input graphs\n";
        return false;
    }
    int useParallelSearch = 0;
    if (isPrimaryProcess())
        useParallelSearch = parallelSearchEligible(graph) ? 1 : 0;
    MPI_Bcast(&useParallelSearch, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (useParallelSearch != 0) return runParallelSearch(graph, output);
#elif defined(ASSEMBLYCPP_USE_OPENMP)
    if (parallelSearchEligible(graph)) return runParallelSearch(graph, output);
#endif

#if defined(ASSEMBLYCPP_USE_MPI)
    int succeeded = 1;
    if (isPrimaryProcess())
    {
        try
        {
            succeeded = improvedBnB(graph, output) ? 1 : 0;
        }
        catch (const std::exception &exception)
        {
            cerr << "error: search failed on the root MPI rank: "
                 << exception.what() << '\n';
            succeeded = 0;
        }
    }
    MPI_Bcast(&succeeded, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return succeeded != 0;
#else
    return improvedBnB(graph, output);
#endif
}

/**
 * @brief Read a molfile or native graph and calculate its assembly index.
 *
 * Molfile paths may include or omit the .mol extension. Native graph paths are
 * used exactly as provided.
 *
 * @param input Input path supplied on the command line.
 * @return true if the input was read and the calculation output was written.
 */
bool assemblyCalculator(const string &input)
{
    searchCancellationFlag.store(false);
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    resetParallelSearchTelemetry();
    resetSearchTelemetry();
#endif
    const bool explicitMolfile = hasMolfileExtension(input);
    const string outputBase =
        explicitMolfile ? input.substr(0, input.size() - 4) : input;
    molGraph mol_graph;
    string inputError;
#if defined(ASSEMBLYCPP_USE_MPI)
    const bool configuredVerbose = verbose;
    if (!isPrimaryProcess()) verbose = false;
#endif
    int inputLoaded = loadMoleculeInput(input, mol_graph, inputError) ? 1 : 0;
#if defined(ASSEMBLYCPP_USE_MPI)
    verbose = configuredVerbose;
    int allInputsLoaded = inputLoaded;
    MPI_Allreduce(
        &inputLoaded,
        &allInputsLoaded,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    inputLoaded = allInputsLoaded;
#endif
    if (!inputLoaded)
    {
        if (isPrimaryProcess())
        {
            if (inputError.empty())
                inputError = "another MPI rank could not load the input";
            cerr << "error: " << inputError << '\n';
        }
        return false;
    }

    const string outputName = outputBase + "Out";
    ofstream outputFile;
    int outputReady = 1;
    if (isPrimaryProcess())
    {
        outputFile.open(outputName);
        if (!outputFile.is_open()) outputReady = 0;
    }
#if defined(ASSEMBLYCPP_USE_MPI)
    MPI_Bcast(&outputReady, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (!outputReady)
    {
        if (isPrimaryProcess())
            cerr << "error: could not open output file '" << outputName << "'\n";
        return false;
    }

    moleculeName = outputBase + "Pathway";
    if (isPrimaryProcess())
        outputFile << outputBase << " has assembly index: ";
    // improvedBnB propagates recoverPathway2's requested-output status.
    bool calculationSucceeded = false;
    try
    {
        calculationSucceeded = runConfiguredSearch(mol_graph, outputFile);
    }
    catch (const std::exception &exception)
    {
        if (isPrimaryProcess())
            cerr << "error: calculation failed for '" << input << "': "
                 << exception.what() << '\n';
        return false;
    }
    bool outputsSucceeded = calculationSucceeded;

    if (
        isPrimaryProcess() && writeIntermediateMAs &&
        !writeoutIntermediateMAs(outputBase + "IntermediateMAs")
    )
    {
        outputsSucceeded = false;
    }

    if (isPrimaryProcess() && receivedUserInterrupt())
    {
        cout << "status: interrupted by user\n";
        outputFile << "status: interrupted by user\n";
    }
    if (isPrimaryProcess())
        outputFile << "time elapsed: " << elapsedClockTicks() << '\n';

    if (isPrimaryProcess()) outputFile.close();
    if (isPrimaryProcess() && !outputFile)
    {
        cerr << "error: could not write output file '" << outputName << "'\n";
        outputsSucceeded = false;
    }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (
        isPrimaryProcess() && searchTelemetryEnabled &&
        !writeSearchTelemetry(outputBase + "Telemetry.json")
    )
    {
        outputsSucceeded = false;
    }
#endif
    return outputsSucceeded;
}
#endif

#ifdef ASSEMBLYCPP_LIBRARY_BUILD

namespace
{

class LibraryOptionScope
{
    int previousEnumerationLimit = ENUM_MAX;
    unsigned long long previousRuntimeTicks = runTimeMax;
    bool previousPathway = isPathway;
    bool previousRemoveHydrogens = removeHydrogens;
    bool previousCompensateDisjoint = disjointCompensation;
    bool previousVerbose = verbose;
    bool previousMemoryReport = memTest;
    bool previousIntermediateMAs = writeIntermediateMAs;

public:
    explicit LibraryOptionScope(const assemblycpp::CalculationOptions &options)
    {
        ENUM_MAX = options.enumerationLimit;
        runTimeMax = options.runtimeTicks;
        isPathway = false;
        removeHydrogens = options.removeHydrogens;
        disjointCompensation = options.compensateDisjoint;
        verbose = options.verbose;
        memTest = false;
        writeIntermediateMAs = false;
    }

    ~LibraryOptionScope()
    {
        ENUM_MAX = previousEnumerationLimit;
        runTimeMax = previousRuntimeTicks;
        isPathway = previousPathway;
        removeHydrogens = previousRemoveHydrogens;
        disjointCompensation = previousCompensateDisjoint;
        verbose = previousVerbose;
        memTest = previousMemoryReport;
        writeIntermediateMAs = previousIntermediateMAs;
    }
};

void prepareInProcessCalculation()
{
    // A runtime budget uses the shared stop flag. Clear it before each item so
    // a limited calculation cannot stop the remainder of an in-process batch.
#ifdef _WIN32
    interruptFlag.store(false);
#else
    interruptFlag = 0;
#endif
    searchCancellationFlag.store(false);
}

assemblycpp::CalculationResult calculateLoadedMolecule(
    molGraph &graph,
    const string &input
)
{
    assemblycpp::CalculationResult result;
    result.input = input;
    prepareInProcessCalculation();
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    resetSearchTelemetry();
#endif
    // An unopened stream is a portable no-file sink for the legacy internal
    // search writer. The public API returns the same value directly below.
    ofstream discardedOutput;
    result.succeeded = improvedBnB(graph, discardedOutput);
    result.assemblyIndex = lastCalculatedAssemblyIndex;
    result.clockTicks = elapsedClockTicks();
    result.runtimeLimitReached = runtimeLimitReached;
    result.enumerationLimitReached = enumerationLimitReached;
    if (!result.succeeded)
        result.error = "calculation did not produce a result";
    return result;
}

bool validLibraryOptions(
    const assemblycpp::CalculationOptions &options,
    string &error
)
{
    if (options.enumerationLimit >= 1) return true;
    error = "enumerationLimit must be at least 1";
    return false;
}

} // namespace

} // namespace assemblycpp::detail
using namespace assemblycpp::detail;

assemblycpp::CalculationResult assemblycpp::calculateMolfile(
    std::istream &molfile,
    const CalculationOptions &options
)
{
    CalculationResult result;
    result.input = "<stream>";
    if (!validLibraryOptions(options, result.error)) return result;

    LibraryOptionScope optionScope(options);
    try
    {
        molGraph graph;
        molfileParser(molfile, graph);
        return calculateLoadedMolecule(graph, result.input);
    }
    catch (const std::exception &exception)
    {
        result.error = exception.what();
        return result;
    }
}

assemblycpp::CalculationResult assemblycpp::calculate(
    const std::string &input,
    const CalculationOptions &options
)
{
    CalculationResult result;
    result.input = input;
    if (!validLibraryOptions(options, result.error)) return result;

    LibraryOptionScope optionScope(options);
    try
    {
        molGraph graph;
        if (!loadMoleculeInput(input, graph, result.error)) return result;
        return calculateLoadedMolecule(graph, input);
    }
    catch (const std::exception &exception)
    {
        result.error = exception.what();
        return result;
    }
}

std::vector<assemblycpp::CalculationResult> assemblycpp::calculateBatch(
    const std::vector<std::string> &inputs,
    const CalculationOptions &options
)
{
    std::vector<CalculationResult> results;
    results.reserve(inputs.size());
    for (const string &input : inputs)
    {
        results.push_back(calculate(input, options));
    }
    return results;
}

#endif

#ifndef ASSEMBLYCPP_LIBRARY_BUILD
/** Write Linux VmPeak memory usage to a file. */
bool maxMemoryUsage(const string& outputFilename)
{
    const string statusFilename = "/proc/self/status";
    ifstream status_file(statusFilename);
    if (!status_file.is_open())
    {
        cerr << "error: could not open memory status file '"
             << statusFilename << "'\n";
        return false;
    }

    string line, peakMemory;

    while (getline(status_file, line))
    {
        if (line.rfind("VmPeak:", 0) == 0)
        {
            peakMemory = line;
            break;
        }
    }

    if (status_file.bad())
    {
        cerr << "error: could not read memory status file '"
             << statusFilename << "'\n";
        return false;
    }
    if (peakMemory.empty())
    {
        cerr << "error: memory status file '" << statusFilename
             << "' does not contain 'VmPeak:'\n";
        return false;
    }

    ofstream outFile(outputFilename);
    if (!outFile.is_open())
    {
        cerr << "error: could not open output file '" << outputFilename << "'\n";
        return false;
    }

    outFile << peakMemory << '\n';
    outFile.close();
    if (!outFile)
    {
        cerr << "error: could not write output file '" << outputFilename << "'\n";
        return false;
    }
    return true;
}
#endif

#ifndef ASSEMBLYCPP_NO_MAIN
#if defined(ASSEMBLYCPP_USE_MPI)
class AssemblyCppMpiSession
{
    bool initialized = false;
    bool usable = false;

public:
    AssemblyCppMpiSession(int &argc, char **&argv)
    {
        int provided = MPI_THREAD_SINGLE;
#if defined(ASSEMBLYCPP_USE_OPENMP)
        constexpr int required = MPI_THREAD_FUNNELED;
#else
        constexpr int required = MPI_THREAD_SINGLE;
#endif
        if (MPI_Init_thread(&argc, &argv, required, &provided) != MPI_SUCCESS)
            return;
        initialized = true;
        MPI_Comm_rank(MPI_COMM_WORLD, &assemblyCppMpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &assemblyCppMpiSize);
        usable = provided >= required;
    }

    ~AssemblyCppMpiSession()
    {
        if (initialized) MPI_Finalize();
    }

    explicit operator bool() const noexcept {return usable;}
};
#endif

int main(int argc, char** argv)
{
    ios::sync_with_stdio(false);
#if defined(ASSEMBLYCPP_USE_MPI)
    AssemblyCppMpiSession mpiSession(argc, argv);
    if (!mpiSession)
    {
        if (assemblyCppMpiRank == 0)
            cerr << "error: MPI could not provide the required thread level\n";
        return 1;
    }
#endif
    #ifdef _WIN32
        static_cast<void>(SetConsoleCtrlHandler(CtrlHandler, TRUE));
    #else
        signal(SIGINT, signalHandler);
    #endif
    CommandLineArguments arguments;
    string argumentError;
    int argumentsValid = 1;
    try
    {
        arguments = parseCommandLine(argc, argv);
    }
    catch (const std::invalid_argument& error)
    {
        argumentError = error.what();
        argumentsValid = 0;
    }

#if defined(ASSEMBLYCPP_USE_MPI)
    int allArgumentsValid = argumentsValid;
    MPI_Allreduce(
        &argumentsValid,
        &allArgumentsValid,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    if (allArgumentsValid == 0)
    {
        if (isPrimaryProcess())
        {
            if (argumentError.empty())
                argumentError = "command line was invalid on another MPI rank";
            cerr << "error: " << argumentError << "\n"
                 << "Usage: AssemblyCpp INPUT [OPTIONS]\n"
                 << "Try 'AssemblyCpp --help' for more information.\n";
        }
        return 2;
    }
    if (!mpiCommandLineOptionsAgree(arguments))
    {
        if (isPrimaryProcess())
            cerr << "error: MPI ranks received different command-line options\n";
        return 2;
    }
#else
    if (argumentsValid == 0)
    {
        cerr << "error: " << argumentError << "\n"
             << "Usage: AssemblyCpp INPUT [OPTIONS]\n"
             << "Try 'AssemblyCpp --help' for more information.\n";
        return 2;
    }
#endif

    if (arguments.showHelp)
    {
        if (isPrimaryProcess()) help();
        return 0;
    }

    bool succeeded = assemblyCalculator(arguments.input);

    // All search and pathway output is complete; ignore new interrupts while
    // final process-level reports and stream buffers are finished.
    disableInterruptHandler();

    #ifdef __linux__
        if (isPrimaryProcess() && succeeded && memTest)
            succeeded = maxMemoryUsage("memUsage");
    #endif

    if (!succeeded) return 1;
    if (receivedUserInterrupt()) return 130;
    return 0;
}
#endif
