#pragma once

/** A root occurrence stored without an EdgeMask from the producer thread. */
struct rootOccurrenceDescriptor
{
    std::size_t wordOffset = 0;
    std::int32_t fragment = 0;
};

/** One independently searchable root matching, expressed only as indices. */
struct rootJobDescriptor
{
    std::size_t firstOccurrence = 0;
    std::size_t secondOccurrence = 0;
    std::int32_t canonicalId = unknownCanonicalId;
    std::uint32_t duplicateSize = 0;
};

static_assert(std::is_trivially_copyable_v<rootOccurrenceDescriptor>);
static_assert(std::is_trivially_copyable_v<rootJobDescriptor>);

/**
 * Mask-free canonical state captured after the one root enumeration.
 *
 * Workers copy this seed into their local canonical caches. The tree interner
 * must travel with graphHashes because tree graph keys contain interner IDs.
 */
struct canonicalisationSeed
{
    decltype(graphHashMap) graphHashes;
    decltype(treeCanonAtomInterner) atomInterner;
    decltype(treeCanonLeafInterner) leafInterner;
    decltype(treeCanonInterner) treeInterner;
};

/**
 * Problem input shared immutably by every worker after construction.
 *
 * In particular, this type deliberately contains no EdgeMask, assemblyState,
 * or assemblyFragment. Wide masks belong to a thread-local arena and must be
 * reconstructed from occurrenceWords by the worker that will destroy them.
 * The two process-local caches are the only mutable members workers touch;
 * both provide their own fine-grained synchronization.
 */
struct SearchContext
{
    molGraph processedMolecule;
    vector<edgeL> universeEdges;
    vector<dagLevel> dag;
    vector<rootOccurrenceDescriptor> rootOccurrences;
    vector<std::uint64_t> occurrenceWords;
    vector<rootJobDescriptor> rootJobs;
    vector<int> homogeneousPathEdgePositions;
    canonicalisationSeed canonicalSeed;
    std::unique_ptr<sharedCanonicalIdRegistry> canonicalRegistry;
    std::unique_ptr<sharedAssemblyTranspositionTable> sharedStates;
    clock_t startedAt = 0;
    unsigned int bondCount = 0;
    int componentCount = 1;
    int rootAssemblyIndex = -1;
    bool enumerationLimit = false;

    SearchContext() = default;
    SearchContext(const SearchContext &) = delete;
    SearchContext &operator=(const SearchContext &) = delete;
    SearchContext(SearchContext &&) = default;
    SearchContext &operator=(SearchContext &&) = default;

    [[nodiscard]] std::size_t maskWordCount() const noexcept
    {
        return universeEdges.size() / EdgeMask::wordBits +
            (universeEdges.size() % EdgeMask::wordBits != 0);
    }
};

struct assemblyPathWitness
{
    vector<assemblyPathStep> current;
    vector<assemblyPathStep> best;
};

inline ASSEMBLYCPP_SEARCH_LOCAL sharedAssemblyTranspositionTable
    *sharedAssemblyStates = nullptr;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t sharedAssemblyWorkerIndex = 0;

/** Mutable caches and fragmentation scratch owned by exactly one worker. */
struct assemblySearchStorage
{
    assemblyTranspositionTable states;
    assemblyPathWitness *pathway = nullptr;
    duplicateClassIndexWorkspace duplicateClassIndex;
    vi candidateKey;

    explicit assemblySearchStorage(assemblyPathWitness *_pathway = nullptr):
        states(1024),
        pathway(_pathway)
    {
        // One search-wide scratch key is safe because every table operation
        // finishes before the synchronous recursive call can reuse it.
        candidateKey.reserve(univEdgeList.size() + 1);
        if (pathway != nullptr)
        {
            pathway->current.reserve(univEdgeList.size());
            pathway->best.reserve(univEdgeList.size());
        }
    }

    /** L1-first lookup for searches configured with a process-shared L2. */
    [[gnu::noinline]] assemblyTranspositionTable::result considerShared(
        std::span<const int> key,
        int sumDupBonds
    )
    {
        const assemblyTranspositionTable::result localResult =
            states.consider(key, sumDupBonds);
        if (localResult == assemblyTranspositionTable::result::dominated)
            return localResult;

        const sharedAssemblyTranspositionTable::consideration sharedResult =
            sharedAssemblyStates->considerWithBestForWorker(
                key,
                sumDupBonds,
                sharedAssemblyWorkerIndex
            );
        if (
            sharedResult.outcome ==
                assemblyTranspositionTable::result::dominated &&
            sharedResult.bestSumDupBonds > sumDupBonds
        )
        {
            // Promote the observed L2 score so subsequent local visits can
            // prune without acquiring the shard again. A later genuinely
            // better score can still improve L1 and be checked in L2.
            static_cast<void>(states.consider(
                key,
                sharedResult.bestSumDupBonds
            ));
        }
        return sharedResult.outcome;
    }
};

/**
 * Per-worker mutable search state. Construct only after configuring that
 * worker's EdgeMask and AtomMask domains from the shared SearchContext.
 */
struct WorkerContext
{
    ufdsMaskWorkspace fragmentation;
    assemblySearchStorage search;
    assemblyState root;
    assemblyState candidate;
    int assemblyIndex;

    explicit WorkerContext(const SearchContext &context):
        fragmentation(
            context.processedMolecule.mg.size(),
            context.universeEdges.size()
        ),
        assemblyIndex(context.rootAssemblyIndex)
    {
        fragmentation.homogeneousPathEdgePositions =
            context.homogeneousPathEdgePositions;

        EdgeMask rootMask;
        rootMask.set();
        root.appendFragment(
            rootMask,
            static_cast<int>(context.universeEdges.size()),
            unknownCanonicalId,
            false
        );
        root.assemblyHashCalculator(search.candidateKey);
        static_cast<void>(search.states.consider(search.candidateKey, 0));
        candidate.reserveFragments(3);
    }

    WorkerContext(const WorkerContext &) = delete;
    WorkerContext &operator=(const WorkerContext &) = delete;
    WorkerContext(WorkerContext &&) = delete;
    WorkerContext &operator=(WorkerContext &&) = delete;
};

/** One fragment of a depth-two task, stored without an owning EdgeMask. */
struct parallelTaskFragmentDescriptor
{
    std::size_t wordOffset = 0;
    std::uint32_t edgeCount = 0;
    bool connected = false;
};

/**
 * A depth-two assembly state that may safely move between worker threads.
 *
 * Canonical IDs are deliberately omitted: after the shared root seed they are
 * process-global only when L2 reuse is enabled, and worker-local otherwise.
 * The receiving worker reconstructs the masks and canonicalises them through
 * the active mode before continuing the search.
 */
struct parallelDepthTwoTaskDescriptor
{
    std::vector<parallelTaskFragmentDescriptor> fragments;
    std::vector<std::uint64_t> fragmentWords;
    int sumDupBonds = 0;
    int lowerBoundAssemblyIndex = 0;
};

/**
 * Rank-local adaptive work controller.
 *
 * Root jobs retain their stable modulo MPI partition and are claimed in
 * leases. When the unclaimed root jobs plus active root work fall below eight
 * per local worker, root searches may donate immediate children until roughly
 * sixteen jobs per worker are ready. The hard policy ceiling is thirty-two
 * queued jobs per worker.
 */
class ParallelTaskScheduler
{
public:
    enum class WorkAvailability
    {
        depthTwo,
        complete,
        wait
    };

    static constexpr std::size_t minimumTasksPerWorker =
        parallelMinimumQueuedTasksPerWorker;
    static constexpr std::size_t targetTasksPerWorker =
        parallelTargetQueuedTasksPerWorker;
    static constexpr std::size_t maximumTasksPerWorker =
        parallelMaximumQueuedTasksPerWorker;
    static constexpr unsigned int maximumTaskDepth =
        parallelMaximumTaskDepth;

    ParallelTaskScheduler(
        std::size_t rootJobCount,
        std::size_t rankPartitionIndex,
        std::size_t rankPartitionCount,
        std::size_t localWorkerCount,
        std::size_t rootLeaseSize,
        bool useAdaptiveRootLeases
    ):
        globalRootJobCount(rootJobCount),
        rankIndex(rankPartitionIndex),
        rankCount(std::max<std::size_t>(1, rankPartitionCount)),
        leaseSize(std::max<std::size_t>(1, rootLeaseSize)),
        rankRootJobCount(
            rootJobCount <= rankIndex
                ? 0
                : 1 + (rootJobCount - 1 - rankIndex) / rankCount
        ),
        workerCount(std::max<std::size_t>(1, localWorkerCount)),
        rootLeaseCursor(0),
        rootJobsOutstanding(rankRootJobCount),
        lowWatermark(saturatedProduct(
            workerCount,
            minimumTasksPerWorker
        )),
        targetWatermark(saturatedProduct(
            workerCount,
            targetTasksPerWorker
        )),
        maximumQueuedTasks(saturatedProduct(
            workerCount,
            maximumTasksPerWorker
        )),
        adaptiveRootLeases(useAdaptiveRootLeases),
        depthTwoEnabled(localWorkerCount > 1),
        proactiveTailRefillEnabled(
            depthTwoEnabled &&
            rankRootJobCount >= saturatedProduct(
                targetWatermark,
                parallelPromisingFrontierLeaseSize
            )
        ),
        idleWorkerTrigger(
            localWorkerCount <= 1
                ? 1
                : (localWorkerCount + 1) / 2
        )
    {
        const std::size_t sparseThreshold = std::max<std::size_t>(
            1,
            lowWatermark / 2
        );
        if (depthTwoEnabled && rankRootJobCount < sparseThreshold)
        {
            refilling = true;
            refillRequested.store(true, std::memory_order_relaxed);
        }
    }

    [[nodiscard]] std::size_t warmStartRootJob() const noexcept
    {
        return globalRootJobCount == 0
            ? std::numeric_limits<std::size_t>::max()
            : 0;
    }

    /** Claim one stable range of rank-partition ordinals. */
    bool claimRootLease(std::size_t &begin, std::size_t &end)
    {
        begin = rootLeaseCursor.load(std::memory_order_relaxed);
        do
        {
            if (begin >= rankRootJobCount) return false;
            const std::size_t claimSize = rootLeaseSize(
                begin,
                rankRootJobCount - begin
            );
            end = claimSize > rankRootJobCount - begin
                ? rankRootJobCount
                : begin + claimSize;
        }
        while (!rootLeaseCursor.compare_exchange_weak(
            begin,
            end,
            std::memory_order_relaxed,
            std::memory_order_relaxed
        ));
        // On substantial frontiers, arm depth-two donation just before the
        // last root tranche so its first child can be handed off. Compact
        // frontiers wait for demonstrated idleness and avoid transfer cost.
        if (
            proactiveTailRefillEnabled &&
            rankRootJobCount - end < lowWatermark &&
            !rootFrontierRefillActivated.load(std::memory_order_relaxed) &&
            !rootFrontierRefillActivated.exchange(
                true,
                std::memory_order_acq_rel
            )
        )
        {
            std::lock_guard<std::mutex> lock(mutex);
            if (!cancelled.load(std::memory_order_relaxed))
            {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                ++searchProactiveTailRefills;
#endif
                refilling = true;
                refillRequested.store(true, std::memory_order_relaxed);
                condition.notify_all();
            }
        }
        return true;
    }

    [[nodiscard]] bool depthTwoRefillRequested() const noexcept
    {
        return refillRequested.load(std::memory_order_relaxed);
    }

    [[nodiscard]] bool depthTwoTasksEnabled() const noexcept
    {
        return depthTwoEnabled;
    }

    /**
     * Queue an immediate child of a root job when the ready frontier is low.
     */
    bool tryEnqueueDepthTwo(
        const assemblyState &state,
        int lowerBoundAssemblyIndex
    )
    {
        if (!depthTwoRefillRequested()) return false;
        {
            std::lock_guard<std::mutex> lock(mutex);
            if (cancelled.load(std::memory_order_relaxed)) return false;

            const std::size_t queued = queuedTaskEstimateLocked();
            if (!refilling)
            {
                refillRequested.store(false, std::memory_order_relaxed);
                return false;
            }
            if (queued >= targetWatermark || queued >= maximumQueuedTasks)
            {
                refilling = false;
                refillRequested.store(false, std::memory_order_relaxed);
                return false;
            }
            ++depthTwoReservations;
            updateRefillStateLocked();
        }

        parallelDepthTwoTaskDescriptor task;
        try
        {
            task.sumDupBonds = state.sumDupBonds;
            task.lowerBoundAssemblyIndex = lowerBoundAssemblyIndex;
            task.fragments.reserve(state.fragments.size());
            const std::size_t wordCount = EdgeMask::activeWordCount();
            if (
                wordCount != 0 &&
                state.fragments.size() >
                    std::numeric_limits<std::size_t>::max() / wordCount
            ) throw std::length_error("depth-two task masks exceed capacity");
            task.fragmentWords.reserve(state.fragments.size() * wordCount);
            for (const assemblyFragment &fragment : state.fragments)
            {
                task.fragments.push_back({
                    task.fragmentWords.size(),
                    fragment.edgeCount,
                    fragment.connected != 0
                });
                for (std::size_t word = 0; word < wordCount; ++word)
                    task.fragmentWords.push_back(fragment.mask.activeWord(word));
            }
        }
        catch (...)
        {
            std::lock_guard<std::mutex> lock(mutex);
            --depthTwoReservations;
            updateRefillStateLocked();
            condition.notify_all();
            throw;
        }

        {
            std::lock_guard<std::mutex> lock(mutex);
            --depthTwoReservations;
            if (cancelled.load(std::memory_order_relaxed))
            {
                updateRefillStateLocked();
                condition.notify_all();
                return false;
            }
            if (!depthTwoQueue.has_value()) depthTwoQueue.emplace();
            depthTwoQueue->push_back(std::move(task));
            ++depthTwoJobsOutstanding;
            updateRefillStateLocked();
            condition.notify_one();
        }
        return true;
    }

    WorkAvailability nextWork(parallelDepthTwoTaskDescriptor &task)
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (depthTwoQueue.has_value() && !depthTwoQueue->empty())
        {
            task = std::move(depthTwoQueue->front());
            depthTwoQueue->pop_front();
            updateRefillStateLocked();
            return WorkAvailability::depthTwo;
        }
        if (
            rootJobsOutstanding.load(std::memory_order_acquire) == 0 &&
            depthTwoJobsOutstanding == 0
        ) return WorkAvailability::complete;
        return WorkAvailability::wait;
    }

    void completeRootLease(std::size_t completedRootJobs)
    {
        if (!depthTwoEnabled) return;
        if (completedRootJobs != 0)
        {
            const std::size_t previous = rootJobsOutstanding.fetch_sub(
                completedRootJobs,
                std::memory_order_acq_rel
            );
            if (previous < completedRootJobs)
            {
                rootJobsOutstanding.fetch_add(
                    completedRootJobs,
                    std::memory_order_relaxed
                );
                throw std::logic_error(
                    "adaptive scheduler completed extra root jobs"
                );
            }
            if (previous == completedRootJobs) condition.notify_all();
        }
    }

    void completeDepthTwoJob()
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (depthTwoJobsOutstanding == 0)
        {
            throw std::logic_error(
                "adaptive scheduler completed an extra depth-two job"
            );
        }
        --depthTwoJobsOutstanding;
        condition.notify_all();
    }

    void waitForWork()
    {
        std::unique_lock<std::mutex> lock(mutex);
        ++idleWorkers;
        const bool workReady = condition.wait_for(
            lock,
            std::chrono::milliseconds(1),
            [&]
            {
                return
                    cancelled.load(std::memory_order_relaxed) ||
                    (depthTwoQueue.has_value() && !depthTwoQueue->empty()) ||
                    rootLeaseCursor.load(std::memory_order_relaxed) <
                        rankRootJobCount ||
                    (rootJobsOutstanding.load(std::memory_order_acquire) == 0 &&
                        depthTwoJobsOutstanding == 0);
            }
        );
        if (
            !workReady &&
            depthTwoEnabled &&
            idleWorkers >= idleWorkerTrigger &&
            !cancelled.load(std::memory_order_relaxed) &&
            rootJobsOutstanding.load(std::memory_order_acquire) != 0 &&
            queuedTaskEstimateLocked() < lowWatermark
        )
        {
            refilling = true;
            refillRequested.store(true, std::memory_order_relaxed);
        }
        --idleWorkers;
    }

    void cancel() noexcept
    {
        cancelled.store(true, std::memory_order_release);
        condition.notify_all();
    }

private:
    static std::size_t saturatedProduct(
        std::size_t left,
        std::size_t right
    ) noexcept
    {
        return left > std::numeric_limits<std::size_t>::max() / right
            ? std::numeric_limits<std::size_t>::max()
            : left * right;
    }

    [[nodiscard]] std::size_t rootLeaseSize(
        std::size_t claimedRootJobs,
        std::size_t remainingRootJobs
    ) const noexcept
    {
        if (!adaptiveRootLeases) return leaseSize;
        // A small initial lease is already fine grained; shrinking it further
        // only adds claim traffic on compact searches.
        if (leaseSize <= parallelPromisingFrontierLeaseSize) return leaseSize;
        if (claimedRootJobs < lowWatermark)
        {
            return std::min(
                std::min(leaseSize, parallelPromisingFrontierLeaseSize),
                remainingRootJobs
            );
        }
        if (remainingRootJobs <= targetWatermark) return 1;
        const std::size_t guided =
            1 + (remainingRootJobs - 1) / targetWatermark;
        return std::min(leaseSize, guided);
    }

    [[nodiscard]] std::size_t queuedTaskEstimateLocked() const noexcept
    {
        const std::size_t claimed = std::min(
            rootLeaseCursor.load(std::memory_order_relaxed),
            rankRootJobCount
        );
        const std::size_t remainingRootJobs = rankRootJobCount - claimed;
        const std::size_t outstanding = rootJobsOutstanding.load(
            std::memory_order_relaxed
        );
        const std::size_t claimedOutstanding = outstanding > remainingRootJobs
            ? outstanding - remainingRootJobs
            : 0;
        const std::size_t activeEstimate = std::min(
            claimedOutstanding,
            workerCount
        );
        if (activeEstimate >
            std::numeric_limits<std::size_t>::max() - remainingRootJobs)
        {
            return std::numeric_limits<std::size_t>::max();
        }
        const std::size_t roots = remainingRootJobs + activeEstimate;
        const std::size_t queuedDepthTwoTasks = depthTwoQueue.has_value()
            ? depthTwoQueue->size()
            : 0;
        if (
            queuedDepthTwoTasks + depthTwoReservations < queuedDepthTwoTasks ||
            queuedDepthTwoTasks + depthTwoReservations >
                std::numeric_limits<std::size_t>::max() - roots
        ) return std::numeric_limits<std::size_t>::max();
        return roots + queuedDepthTwoTasks + depthTwoReservations;
    }

    void updateRefillStateLocked() noexcept
    {
        const std::size_t queued = queuedTaskEstimateLocked();
        if (refilling && queued >= targetWatermark) refilling = false;
        refillRequested.store(refilling, std::memory_order_relaxed);
    }

    const std::size_t globalRootJobCount;
    const std::size_t rankIndex;
    const std::size_t rankCount;
    const std::size_t leaseSize;
    const std::size_t rankRootJobCount;
    const std::size_t workerCount;
    std::atomic<std::size_t> rootLeaseCursor;
    mutable std::mutex mutex;
    std::condition_variable condition;
    std::optional<std::deque<parallelDepthTwoTaskDescriptor>> depthTwoQueue;
    std::size_t depthTwoReservations = 0;
    std::atomic<std::size_t> rootJobsOutstanding;
    std::size_t depthTwoJobsOutstanding = 0;
    const std::size_t lowWatermark;
    const std::size_t targetWatermark;
    const std::size_t maximumQueuedTasks;
    const bool adaptiveRootLeases;
    const bool depthTwoEnabled;
    const bool proactiveTailRefillEnabled;
    const std::size_t idleWorkerTrigger;
    std::size_t idleWorkers = 0;
    bool refilling = false;
    std::atomic_bool rootFrontierRefillActivated{false};
    std::atomic_bool refillRequested{false};
    std::atomic_bool cancelled{false};
};

ASSEMBLYCPP_SEARCH_LOCAL ParallelTaskScheduler *parallelTaskScheduler = nullptr;
ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchDepthTwoTasksSpawned = 0;
ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchDepthTwoTasksExecuted = 0;
ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchWarmStartBranches = 0;
