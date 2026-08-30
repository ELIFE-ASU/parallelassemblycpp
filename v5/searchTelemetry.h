#pragma once

#ifndef ASSEMBLY_ENABLE_TELEMETRY

inline constexpr bool searchTelemetryCompiled = false;
inline bool searchTelemetryEnabled = false;

#else

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef __linux__
    #include <cerrno>
    #include <cstdlib>
    #include <cstring>
    #include <fcntl.h>
    #include <unistd.h>
#endif

enum class SearchTelemetryPhase : size_t
{
    inputSetup,
    initialEnumeration,
    dagConversion,
    assemblySearch,
    output,
    count
};

struct SearchTelemetryCounters
{
    uint64_t retainedMaskAttempts = 0;
    uint64_t retainedMasks = 0;
    uint64_t duplicateMaskAttempts = 0;
    uint64_t rejectedMasks = 0;
    uint64_t matchingVisits = 0;

    uint64_t canonicalisationCalls = 0;
    uint64_t canonicalisationMaskCacheHits = 0;
    uint64_t canonicalisationMaskCacheMisses = 0;
    uint64_t canonicalClassInsertions = 0;
    uint64_t canonicalClassReuses = 0;
    uint64_t vf2Calls = 0;
    uint64_t vf2Matches = 0;

    uint64_t residualDecompositionRequests = 0;
    uint64_t residualCacheEligibleRequests = 0;
    uint64_t residualCacheSmallMoleculeBypasses = 0;
    uint64_t residualCacheWideMoleculeBypasses = 0;
    uint64_t residualCacheSmallResidualBypasses = 0;
    uint64_t residualCacheFirstOccurrenceBypasses = 0;
    uint64_t residualCacheRuntimeDisabledBypasses = 0;
    uint64_t residualCacheLookups = 0;
    uint64_t residualCacheHits = 0;
    uint64_t residualCacheMisses = 0;
    uint64_t residualCacheAdmissions = 0;

    uint64_t assemblyCacheLookups = 0;
    uint64_t assemblyCacheHits = 0;
    uint64_t assemblyCacheMisses = 0;
    uint64_t assemblyCachePrunedHits = 0;
    uint64_t assemblyCacheUpdatedHits = 0;

    uint64_t pairBoundCacheLookups = 0;
    uint64_t pairBoundCacheHits = 0;
    uint64_t pairBoundCacheMisses = 0;
};

struct ProcessMemorySnapshot
{
    bool residentAvailable = false;
    bool virtualAvailable = false;
    uint64_t residentKiB = 0;
    uint64_t residentHighWaterKiB = 0;
    uint64_t virtualKiB = 0;
    uint64_t virtualPeakKiB = 0;
};

struct SearchTelemetryPhaseStats
{
    uint64_t clockTicks = 0;
    uint64_t wallNanoseconds = 0;
    uint64_t startResidentKiB = 0;
    uint64_t peakResidentKiB = 0;
    uint64_t endResidentKiB = 0;
    uint64_t startVirtualKiB = 0;
    uint64_t endVirtualKiB = 0;
    size_t activations = 0;
    bool startResidentAvailable = false;
    bool endResidentAvailable = false;
    bool startVirtualAvailable = false;
    bool endVirtualAvailable = false;
    bool exactResidentPeak = false;
};

struct SearchTelemetryState
{
    SearchTelemetryCounters counters;
    std::array<
        SearchTelemetryPhaseStats,
        static_cast<size_t>(SearchTelemetryPhase::count)
    > phases;
    SearchTelemetryPhase currentPhase = SearchTelemetryPhase::count;
    clock_t phaseStart = 0;
    uint64_t phaseWallStartNanoseconds = 0;
    size_t processedAtoms = 0;
    size_t processedEdges = 0;
    size_t activeMaskWords = 0;
    bool residualCacheEligible = false;
    bool collectPhaseMemory = true;
    bool active = false;
    ProcessMemorySnapshot finalMemory;
};

/** One fixed-width worker record suitable for process and MPI reduction. */
struct ParallelSearchWorkerTelemetry
{
    uint64_t mpiRank = 0;
    uint64_t localWorkerIndex = 0;
    uint64_t globalWorkerIndex = 0;
    uint64_t rankPartitionIndex = 0;
    uint64_t rankPartitionCount = 1;
    uint64_t branchCandidates = 0;
    uint64_t branchLeases = 0;
    uint64_t branchAssignments = 0;
    uint64_t depthTwoTasksSpawned = 0;
    uint64_t depthTwoTasksExecuted = 0;
    uint64_t deeperTasksSpawned = 0;
    uint64_t deeperTasksExecuted = 0;
    uint64_t taskStealAttempts = 0;
    uint64_t taskSteals = 0;
    uint64_t localTaskExecutions = 0;
    uint64_t schedulerIdleWaits = 0;
    uint64_t schedulerIdleNanoseconds = 0;
    uint64_t deepRefillActivations = 0;
    uint64_t taskQueueHighWatermark = 0;
    uint64_t maximumTaskDepthExecuted = 0;
    uint64_t proactiveTailRefills = 0;
    uint64_t warmStartBranches = 0;
    uint64_t elapsedNanoseconds = 0;
    uint64_t busyNanoseconds = 0;
    uint64_t processedAtoms = 0;
    uint64_t processedEdges = 0;
    uint64_t activeMaskWords = 0;
    uint64_t residualCacheEligible = 0;
    SearchTelemetryCounters counters;
    std::array<
        uint64_t,
        static_cast<size_t>(SearchTelemetryPhase::count)
    > phaseClockTicks{};
    std::array<
        uint64_t,
        static_cast<size_t>(SearchTelemetryPhase::count)
    > phaseWallNanoseconds{};
    std::array<
        uint64_t,
        static_cast<size_t>(SearchTelemetryPhase::count)
    > phaseActivations{};
};

static_assert(std::is_trivially_copyable_v<ParallelSearchWorkerTelemetry>);

struct ParallelSearchTelemetrySummary
{
    bool enabled = false;
    bool branchScanComplete = false;
    bool branchSchedulerComplete = false;
    bool branchCandidateCountsAgree = false;
    std::string mode;
    std::string aggregationScope;
    uint64_t rankCount = 1;
    uint64_t workerCount = 1;
    uint64_t branchLeaseSize = 1;
    uint64_t branchCandidateCount = 0;
    uint64_t branchLeaseCount = 0;
    uint64_t branchAssignmentCount = 0;
    uint64_t depthTwoTaskSpawnCount = 0;
    uint64_t depthTwoTaskExecutionCount = 0;
    uint64_t deeperTaskSpawnCount = 0;
    uint64_t deeperTaskExecutionCount = 0;
    uint64_t taskStealAttemptCount = 0;
    uint64_t taskStealCount = 0;
    uint64_t localTaskExecutionCount = 0;
    uint64_t schedulerIdleWaitCount = 0;
    uint64_t schedulerIdleNanoseconds = 0;
    uint64_t deepRefillActivationCount = 0;
    uint64_t taskQueueHighWatermark = 0;
    uint64_t maximumTaskDepthExecuted = 0;
    uint64_t proactiveTailRefillCount = 0;
    uint64_t warmStartBranchCount = 0;
    uint64_t elapsedNanoseconds = 0;
    uint64_t workerElapsedNanoseconds = 0;
    uint64_t workerBusyNanoseconds = 0;
    std::vector<uint64_t> localThreadsPerRank;
    SearchTelemetryCounters aggregateCounters;
    std::vector<ParallelSearchWorkerTelemetry> workers;
};

inline constexpr bool searchTelemetryCompiled = true;

inline bool searchTelemetryEnabled = false;
inline ASSEMBLYCPP_SEARCH_LOCAL SearchTelemetryState searchTelemetry;
inline ParallelSearchTelemetrySummary parallelSearchTelemetry;

inline const char* searchTelemetryPhaseName(SearchTelemetryPhase phase)
{
    switch (phase)
    {
        case SearchTelemetryPhase::inputSetup: return "input_setup";
        case SearchTelemetryPhase::initialEnumeration:
            return "initial_enumeration";
        case SearchTelemetryPhase::dagConversion: return "dag_conversion";
        case SearchTelemetryPhase::assemblySearch: return "assembly_search";
        case SearchTelemetryPhase::output: return "output";
        case SearchTelemetryPhase::count: break;
    }
    return "unknown";
}

inline uint64_t telemetryClockDifference(clock_t start, clock_t end)
{
    const clock_t clockError = static_cast<clock_t>(-1);
    if (start == clockError || end == clockError) return 0;
    using UnsignedClock = std::make_unsigned_t<clock_t>;
    return static_cast<uint64_t>(
        static_cast<UnsignedClock>(end) - static_cast<UnsignedClock>(start)
    );
}

inline uint64_t searchTelemetryWallNanoseconds()
{
    using namespace std::chrono;
    return static_cast<uint64_t>(duration_cast<nanoseconds>(
        steady_clock::now().time_since_epoch()
    ).count());
}

inline uint64_t telemetryNanosecondDifference(uint64_t start, uint64_t end)
{
    return end >= start ? end - start : 0;
}

#ifdef __linux__
inline bool parseProcStatusValue(
    const char *buffer,
    const char *end,
    const char *label,
    uint64_t &value
)
{
    const size_t labelLength = std::strlen(label);
    const char *position = buffer;
    while (position < end)
    {
        const char *lineEnd = static_cast<const char *>(
            std::memchr(position, '\n', static_cast<size_t>(end - position))
        );
        if (lineEnd == nullptr) lineEnd = end;
        if (
            static_cast<size_t>(lineEnd - position) >= labelLength &&
            std::memcmp(position, label, labelLength) == 0
        )
        {
            const char *number = position + labelLength;
            while (number < lineEnd && (*number == ' ' || *number == '\t'))
                ++number;
            errno = 0;
            char *parsedEnd = nullptr;
            const unsigned long long parsed = std::strtoull(
                number,
                &parsedEnd,
                10
            );
            if (
                errno != 0 ||
                parsedEnd == number ||
                parsedEnd > lineEnd
            ) return false;
            value = static_cast<uint64_t>(parsed);
            return true;
        }
        position = lineEnd < end ? lineEnd + 1 : end;
    }
    return false;
}

inline ProcessMemorySnapshot readProcessMemorySnapshot()
{
    ProcessMemorySnapshot result;
    const int descriptor = open("/proc/self/status", O_RDONLY | O_CLOEXEC);
    if (descriptor < 0) return result;

    std::array<char, 16384> buffer{};
    size_t used = 0;
    while (used + 1 < buffer.size())
    {
        const ssize_t count = read(
            descriptor,
            buffer.data() + used,
            buffer.size() - used - 1
        );
        if (count == 0) break;
        if (count < 0)
        {
            if (errno == EINTR) continue;
            close(descriptor);
            return result;
        }
        used += static_cast<size_t>(count);
    }
    close(descriptor);
    const char *end = buffer.data() + used;
    result.residentAvailable =
        parseProcStatusValue(
            buffer.data(), end, "VmRSS:", result.residentKiB
        ) &&
        parseProcStatusValue(
            buffer.data(), end, "VmHWM:", result.residentHighWaterKiB
        );
    result.virtualAvailable =
        parseProcStatusValue(
            buffer.data(), end, "VmSize:", result.virtualKiB
        ) &&
        parseProcStatusValue(
            buffer.data(), end, "VmPeak:", result.virtualPeakKiB
        );
    return result;
}

inline bool resetProcessResidentHighWaterMark()
{
    const int descriptor = open("/proc/self/clear_refs", O_WRONLY | O_CLOEXEC);
    if (descriptor < 0) return false;
    constexpr char resetCommand[] = "5\n";
    ssize_t written = 0;
    do
    {
        written = write(descriptor, resetCommand, sizeof(resetCommand) - 1);
    } while (written < 0 && errno == EINTR);
    const bool closed = close(descriptor) == 0;
    return written == static_cast<ssize_t>(sizeof(resetCommand) - 1) && closed;
}
#else
inline ProcessMemorySnapshot readProcessMemorySnapshot()
{
    return {};
}

inline bool resetProcessResidentHighWaterMark()
{
    return false;
}
#endif

inline void finishSearchTelemetryPhase()
{
    if (
        !searchTelemetry.active ||
        searchTelemetry.currentPhase == SearchTelemetryPhase::count
    ) return;

    const clock_t now = clock();
    const uint64_t wallNow = searchTelemetryWallNanoseconds();
    SearchTelemetryPhaseStats &stats = searchTelemetry.phases[
        static_cast<size_t>(searchTelemetry.currentPhase)
    ];
    stats.clockTicks += telemetryClockDifference(
        searchTelemetry.phaseStart,
        now
    );
    stats.wallNanoseconds += telemetryNanosecondDifference(
        searchTelemetry.phaseWallStartNanoseconds,
        wallNow
    );
    if (!searchTelemetry.collectPhaseMemory) return;
    ProcessMemorySnapshot memory = readProcessMemorySnapshot();
    stats.endResidentAvailable = memory.residentAvailable;
    if (memory.residentAvailable)
    {
        stats.endResidentKiB = memory.residentKiB;
        if (stats.exactResidentPeak)
        {
            if (memory.residentHighWaterKiB > stats.peakResidentKiB)
                stats.peakResidentKiB = memory.residentHighWaterKiB;
        }
    }
    else stats.exactResidentPeak = false;
    stats.endVirtualAvailable = memory.virtualAvailable;
    if (memory.virtualAvailable)
    {
        stats.endVirtualKiB = memory.virtualKiB;
    }
}

inline void setSearchTelemetryPhase(SearchTelemetryPhase phase)
{
    if (
        !searchTelemetryCompiled ||
        !searchTelemetryEnabled ||
        !searchTelemetry.active ||
        phase == SearchTelemetryPhase::count
    ) return;
    if (searchTelemetry.currentPhase == phase) return;
    finishSearchTelemetryPhase();

    SearchTelemetryPhaseStats &stats = searchTelemetry.phases[
        static_cast<size_t>(phase)
    ];
    const bool firstActivation = stats.activations == 0;
    ++stats.activations;
    if (searchTelemetry.collectPhaseMemory)
    {
        const bool exactResidentPeak = resetProcessResidentHighWaterMark();
        ProcessMemorySnapshot memory = readProcessMemorySnapshot();
        stats.exactResidentPeak =
            (firstActivation || stats.exactResidentPeak) &&
            exactResidentPeak &&
            memory.residentAvailable;
        stats.endResidentAvailable = memory.residentAvailable;
        if (memory.residentAvailable)
        {
            if (firstActivation)
            {
                stats.startResidentAvailable = true;
                stats.startResidentKiB = memory.residentKiB;
            }
            stats.endResidentKiB = memory.residentKiB;
            if (exactResidentPeak)
            {
                if (memory.residentHighWaterKiB > stats.peakResidentKiB)
                    stats.peakResidentKiB = memory.residentHighWaterKiB;
            }
        }
        stats.endVirtualAvailable = memory.virtualAvailable;
        if (memory.virtualAvailable)
        {
            if (firstActivation)
            {
                stats.startVirtualAvailable = true;
                stats.startVirtualKiB = memory.virtualKiB;
            }
            stats.endVirtualKiB = memory.virtualKiB;
        }
    }
    searchTelemetry.currentPhase = phase;
    searchTelemetry.phaseStart = clock();
    searchTelemetry.phaseWallStartNanoseconds =
        searchTelemetryWallNanoseconds();
}

inline void resetSearchTelemetry(bool collectPhaseMemory = true)
{
    if (!searchTelemetryCompiled || !searchTelemetryEnabled) return;
    searchTelemetry = SearchTelemetryState{};
    searchTelemetry.collectPhaseMemory = collectPhaseMemory;
    searchTelemetry.active = true;
    setSearchTelemetryPhase(SearchTelemetryPhase::inputSetup);
}

inline void configureSearchTelemetryGraph(
    size_t atoms,
    size_t edges,
    size_t activeMaskWords,
    bool residualCacheEligible
)
{
    if (!searchTelemetryCompiled || !searchTelemetryEnabled) return;
    searchTelemetry.processedAtoms = atoms;
    searchTelemetry.processedEdges = edges;
    searchTelemetry.activeMaskWords = activeMaskWords;
    searchTelemetry.residualCacheEligible = residualCacheEligible;
}

inline void finaliseSearchTelemetry()
{
    if (
        !searchTelemetryCompiled ||
        !searchTelemetryEnabled ||
        !searchTelemetry.active
    ) return;
    finishSearchTelemetryPhase();
    if (searchTelemetry.collectPhaseMemory)
        searchTelemetry.finalMemory = readProcessMemorySnapshot();
    searchTelemetry.currentPhase = SearchTelemetryPhase::count;
    searchTelemetry.active = false;
}

inline void addSearchTelemetryCounters(
    SearchTelemetryCounters &destination,
    const SearchTelemetryCounters &source
)
{
    destination.retainedMaskAttempts += source.retainedMaskAttempts;
    destination.retainedMasks += source.retainedMasks;
    destination.duplicateMaskAttempts += source.duplicateMaskAttempts;
    destination.rejectedMasks += source.rejectedMasks;
    destination.matchingVisits += source.matchingVisits;

    destination.canonicalisationCalls += source.canonicalisationCalls;
    destination.canonicalisationMaskCacheHits +=
        source.canonicalisationMaskCacheHits;
    destination.canonicalisationMaskCacheMisses +=
        source.canonicalisationMaskCacheMisses;
    destination.canonicalClassInsertions +=
        source.canonicalClassInsertions;
    destination.canonicalClassReuses += source.canonicalClassReuses;
    destination.vf2Calls += source.vf2Calls;
    destination.vf2Matches += source.vf2Matches;

    destination.residualDecompositionRequests +=
        source.residualDecompositionRequests;
    destination.residualCacheEligibleRequests +=
        source.residualCacheEligibleRequests;
    destination.residualCacheSmallMoleculeBypasses +=
        source.residualCacheSmallMoleculeBypasses;
    destination.residualCacheWideMoleculeBypasses +=
        source.residualCacheWideMoleculeBypasses;
    destination.residualCacheSmallResidualBypasses +=
        source.residualCacheSmallResidualBypasses;
    destination.residualCacheFirstOccurrenceBypasses +=
        source.residualCacheFirstOccurrenceBypasses;
    destination.residualCacheRuntimeDisabledBypasses +=
        source.residualCacheRuntimeDisabledBypasses;
    destination.residualCacheLookups += source.residualCacheLookups;
    destination.residualCacheHits += source.residualCacheHits;
    destination.residualCacheMisses += source.residualCacheMisses;
    destination.residualCacheAdmissions += source.residualCacheAdmissions;

    destination.assemblyCacheLookups += source.assemblyCacheLookups;
    destination.assemblyCacheHits += source.assemblyCacheHits;
    destination.assemblyCacheMisses += source.assemblyCacheMisses;
    destination.assemblyCachePrunedHits += source.assemblyCachePrunedHits;
    destination.assemblyCacheUpdatedHits += source.assemblyCacheUpdatedHits;

    destination.pairBoundCacheLookups += source.pairBoundCacheLookups;
    destination.pairBoundCacheHits += source.pairBoundCacheHits;
    destination.pairBoundCacheMisses += source.pairBoundCacheMisses;
}

inline ParallelSearchWorkerTelemetry captureParallelSearchWorkerTelemetry(
    uint64_t mpiRank,
    uint64_t localWorkerIndex,
    uint64_t globalWorkerIndex,
    uint64_t rankPartitionIndex,
    uint64_t rankPartitionCount,
    uint64_t branchCandidates,
    uint64_t branchLeases,
    uint64_t branchAssignments,
    uint64_t depthTwoTasksSpawned,
    uint64_t depthTwoTasksExecuted,
    uint64_t deeperTasksSpawned,
    uint64_t deeperTasksExecuted,
    uint64_t taskStealAttempts,
    uint64_t taskSteals,
    uint64_t localTaskExecutions,
    uint64_t schedulerIdleWaits,
    uint64_t schedulerIdleNanoseconds,
    uint64_t deepRefillActivations,
    uint64_t taskQueueHighWatermark,
    uint64_t maximumTaskDepthExecuted,
    uint64_t proactiveTailRefills,
    uint64_t warmStartBranches,
    uint64_t elapsedNanoseconds
)
{
    finaliseSearchTelemetry();
    ParallelSearchWorkerTelemetry result;
    result.mpiRank = mpiRank;
    result.localWorkerIndex = localWorkerIndex;
    result.globalWorkerIndex = globalWorkerIndex;
    result.rankPartitionIndex = rankPartitionIndex;
    result.rankPartitionCount = rankPartitionCount;
    result.branchCandidates = branchCandidates;
    result.branchLeases = branchLeases;
    result.branchAssignments = branchAssignments;
    result.depthTwoTasksSpawned = depthTwoTasksSpawned;
    result.depthTwoTasksExecuted = depthTwoTasksExecuted;
    result.deeperTasksSpawned = deeperTasksSpawned;
    result.deeperTasksExecuted = deeperTasksExecuted;
    result.taskStealAttempts = taskStealAttempts;
    result.taskSteals = taskSteals;
    result.localTaskExecutions = localTaskExecutions;
    result.schedulerIdleWaits = schedulerIdleWaits;
    result.schedulerIdleNanoseconds = schedulerIdleNanoseconds;
    result.deepRefillActivations = deepRefillActivations;
    result.taskQueueHighWatermark = taskQueueHighWatermark;
    result.maximumTaskDepthExecuted = maximumTaskDepthExecuted;
    result.proactiveTailRefills = proactiveTailRefills;
    result.warmStartBranches = warmStartBranches;
    result.elapsedNanoseconds = elapsedNanoseconds;
    result.processedAtoms = searchTelemetry.processedAtoms;
    result.processedEdges = searchTelemetry.processedEdges;
    result.activeMaskWords = searchTelemetry.activeMaskWords;
    result.residualCacheEligible =
        searchTelemetry.residualCacheEligible ? 1 : 0;
    result.counters = searchTelemetry.counters;
    for (size_t i = 0; i < result.phaseClockTicks.size(); ++i)
    {
        const SearchTelemetryPhaseStats &phase = searchTelemetry.phases[i];
        result.phaseClockTicks[i] = phase.clockTicks;
        result.phaseWallNanoseconds[i] = phase.wallNanoseconds;
        result.phaseActivations[i] = static_cast<uint64_t>(phase.activations);
    }
    result.busyNanoseconds = result.elapsedNanoseconds - std::min(
        result.schedulerIdleNanoseconds,
        result.elapsedNanoseconds
    );
    return result;
}

inline void resetParallelSearchTelemetry()
{
    parallelSearchTelemetry = ParallelSearchTelemetrySummary{};
}

inline void configureParallelSearchTelemetry(
    std::string mode,
    std::string aggregationScope,
    uint64_t rankCount,
    uint64_t branchLeaseSize,
    uint64_t elapsedNanoseconds,
    bool completedSearch,
    std::vector<ParallelSearchWorkerTelemetry> workers
)
{
    resetParallelSearchTelemetry();
    if (workers.empty()) return;

    std::sort(
        workers.begin(),
        workers.end(),
        [](const auto &left, const auto &right)
        {
            return left.globalWorkerIndex < right.globalWorkerIndex;
        }
    );

    ParallelSearchTelemetrySummary &summary = parallelSearchTelemetry;
    summary.enabled = true;
    summary.mode = std::move(mode);
    summary.aggregationScope = std::move(aggregationScope);
    summary.rankCount = std::max<uint64_t>(1, rankCount);
    summary.workerCount = static_cast<uint64_t>(workers.size());
    summary.branchLeaseSize = branchLeaseSize;
    summary.elapsedNanoseconds = elapsedNanoseconds;
    summary.localThreadsPerRank.assign(summary.rankCount, 0);
    summary.branchSchedulerComplete = branchLeaseSize > 0;
    summary.branchCandidateCountsAgree = true;
    summary.branchCandidateCount = workers.front().branchCandidates;

    std::vector<bool> seenGlobalWorker(summary.workerCount, false);
    std::vector<uint64_t> branchAssignmentsPerRank(summary.rankCount, 0);
    const ParallelSearchWorkerTelemetry &first = workers.front();
    SearchTelemetryState merged;
    merged.collectPhaseMemory = false;
    merged.processedAtoms = static_cast<size_t>(first.processedAtoms);
    merged.processedEdges = static_cast<size_t>(first.processedEdges);
    merged.activeMaskWords = static_cast<size_t>(first.activeMaskWords);
    merged.residualCacheEligible = first.residualCacheEligible != 0;

    for (const ParallelSearchWorkerTelemetry &worker : workers)
    {
        if (worker.mpiRank < summary.rankCount)
        {
            ++summary.localThreadsPerRank[worker.mpiRank];
            branchAssignmentsPerRank[worker.mpiRank] +=
                worker.branchAssignments;
        }
        else summary.branchSchedulerComplete = false;

        if (
            worker.globalWorkerIndex >= summary.workerCount ||
            seenGlobalWorker[worker.globalWorkerIndex]
        ) summary.branchSchedulerComplete = false;
        else seenGlobalWorker[worker.globalWorkerIndex] = true;

        if (
            worker.rankPartitionIndex != worker.mpiRank ||
            worker.rankPartitionCount != summary.rankCount
        ) summary.branchSchedulerComplete = false;

        if (
            (worker.branchAssignments == 0 && worker.branchLeases != 0) ||
            (worker.branchAssignments > 0 && (
                summary.branchLeaseSize == 0 ||
                worker.branchLeases == 0 ||
                worker.branchLeases > worker.branchAssignments ||
                1 + (worker.branchAssignments - 1) / summary.branchLeaseSize >
                    worker.branchLeases
            ))
        ) summary.branchSchedulerComplete = false;

        if (worker.branchCandidates != summary.branchCandidateCount)
            summary.branchCandidateCountsAgree = false;
        if (
            worker.processedAtoms != first.processedAtoms ||
            worker.processedEdges != first.processedEdges ||
            worker.activeMaskWords != first.activeMaskWords ||
            worker.residualCacheEligible != first.residualCacheEligible
        ) summary.branchCandidateCountsAgree = false;

        summary.branchLeaseCount += worker.branchLeases;
        summary.branchAssignmentCount += worker.branchAssignments;
        summary.depthTwoTaskSpawnCount += worker.depthTwoTasksSpawned;
        summary.depthTwoTaskExecutionCount += worker.depthTwoTasksExecuted;
        summary.deeperTaskSpawnCount += worker.deeperTasksSpawned;
        summary.deeperTaskExecutionCount += worker.deeperTasksExecuted;
        summary.taskStealAttemptCount += worker.taskStealAttempts;
        summary.taskStealCount += worker.taskSteals;
        summary.localTaskExecutionCount += worker.localTaskExecutions;
        summary.schedulerIdleWaitCount += worker.schedulerIdleWaits;
        summary.schedulerIdleNanoseconds += worker.schedulerIdleNanoseconds;
        summary.deepRefillActivationCount += worker.deepRefillActivations;
        summary.taskQueueHighWatermark = std::max(
            summary.taskQueueHighWatermark,
            worker.taskQueueHighWatermark
        );
        summary.maximumTaskDepthExecuted = std::max(
            summary.maximumTaskDepthExecuted,
            worker.maximumTaskDepthExecuted
        );
        summary.proactiveTailRefillCount += worker.proactiveTailRefills;
        summary.warmStartBranchCount += worker.warmStartBranches;
        summary.workerElapsedNanoseconds += worker.elapsedNanoseconds;
        summary.workerBusyNanoseconds += worker.busyNanoseconds;
        addSearchTelemetryCounters(summary.aggregateCounters, worker.counters);

        for (size_t i = 0; i < merged.phases.size(); ++i)
        {
            merged.phases[i].wallNanoseconds +=
                worker.phaseWallNanoseconds[i];
            merged.phases[i].activations +=
                static_cast<size_t>(worker.phaseActivations[i]);
        }
    }
    // clock() is process CPU time. Keep one observation for the legacy field;
    // summing overlapping worker observations would overstate CPU time.
    for (size_t i = 0; i < merged.phases.size(); ++i)
        merged.phases[i].clockTicks = first.phaseClockTicks[i];

    if (summary.branchCandidateCountsAgree)
    {
        for (uint64_t rank = 0; rank < summary.rankCount; ++rank)
        {
            const uint64_t expectedAssignments =
                summary.branchCandidateCount <= rank
                ? 0
                : 1 +
                    (summary.branchCandidateCount - 1 - rank) /
                    summary.rankCount;
            if (branchAssignmentsPerRank[rank] != expectedAssignments)
                summary.branchSchedulerComplete = false;
        }
    }
    const uint64_t transferredTaskExecutions =
        summary.depthTwoTaskExecutionCount +
        summary.deeperTaskExecutionCount;
    summary.branchSchedulerComplete =
        completedSearch &&
        summary.branchSchedulerComplete &&
        summary.branchCandidateCountsAgree &&
        summary.branchAssignmentCount == summary.branchCandidateCount &&
        summary.depthTwoTaskSpawnCount ==
            summary.depthTwoTaskExecutionCount &&
        summary.deeperTaskSpawnCount == summary.deeperTaskExecutionCount &&
        summary.taskStealCount <= summary.taskStealAttemptCount &&
        summary.localTaskExecutionCount + summary.taskStealCount ==
            transferredTaskExecutions &&
        (
            (transferredTaskExecutions == 0 &&
                summary.maximumTaskDepthExecuted == 0) ||
            (transferredTaskExecutions != 0 &&
                summary.maximumTaskDepthExecuted >= 2 &&
                summary.maximumTaskDepthExecuted <=
                    parallelMaximumTaskDepth)
        ) &&
        summary.warmStartBranchCount ==
            (summary.branchCandidateCount == 0 ? 0 : summary.rankCount);
    summary.branchScanComplete = summary.branchSchedulerComplete;
    summary.workers = std::move(workers);
    merged.counters = summary.aggregateCounters;
    searchTelemetry = std::move(merged);
}

inline void writeAllSearchTelemetryCounters(
    std::ostream &output,
    const SearchTelemetryCounters &counters,
    const std::string &indent
)
{
    output << "{\n"
           << indent << "  \"retained_mask_attempts\": "
           << counters.retainedMaskAttempts << ",\n"
           << indent << "  \"retained_masks\": "
           << counters.retainedMasks << ",\n"
           << indent << "  \"duplicate_mask_attempts\": "
           << counters.duplicateMaskAttempts << ",\n"
           << indent << "  \"rejected_masks\": "
           << counters.rejectedMasks << ",\n"
           << indent << "  \"matching_visits\": "
           << counters.matchingVisits << ",\n"
           << indent << "  \"canonicalisation_calls\": "
           << counters.canonicalisationCalls << ",\n"
           << indent << "  \"canonicalisation_mask_cache_hits\": "
           << counters.canonicalisationMaskCacheHits << ",\n"
           << indent << "  \"canonicalisation_mask_cache_misses\": "
           << counters.canonicalisationMaskCacheMisses << ",\n"
           << indent << "  \"canonical_class_insertions\": "
           << counters.canonicalClassInsertions << ",\n"
           << indent << "  \"canonical_class_reuses\": "
           << counters.canonicalClassReuses << ",\n"
           << indent << "  \"vf2_calls\": " << counters.vf2Calls << ",\n"
           << indent << "  \"vf2_matches\": " << counters.vf2Matches
           << ",\n"
           << indent << "  \"residual_decomposition_requests\": "
           << counters.residualDecompositionRequests << ",\n"
           << indent << "  \"residual_cache_eligible_requests\": "
           << counters.residualCacheEligibleRequests << ",\n"
           << indent << "  \"residual_cache_small_molecule_bypasses\": "
           << counters.residualCacheSmallMoleculeBypasses << ",\n"
           << indent << "  \"residual_cache_wide_molecule_bypasses\": "
           << counters.residualCacheWideMoleculeBypasses << ",\n"
           << indent << "  \"residual_cache_small_residual_bypasses\": "
           << counters.residualCacheSmallResidualBypasses << ",\n"
           << indent << "  \"residual_cache_first_occurrence_bypasses\": "
           << counters.residualCacheFirstOccurrenceBypasses << ",\n"
           << indent << "  \"residual_cache_runtime_disabled_bypasses\": "
           << counters.residualCacheRuntimeDisabledBypasses << ",\n"
           << indent << "  \"residual_cache_lookups\": "
           << counters.residualCacheLookups << ",\n"
           << indent << "  \"residual_cache_hits\": "
           << counters.residualCacheHits << ",\n"
           << indent << "  \"residual_cache_misses\": "
           << counters.residualCacheMisses << ",\n"
           << indent << "  \"residual_cache_admissions\": "
           << counters.residualCacheAdmissions << ",\n"
           << indent << "  \"assembly_cache_lookups\": "
           << counters.assemblyCacheLookups << ",\n"
           << indent << "  \"assembly_cache_hits\": "
           << counters.assemblyCacheHits << ",\n"
           << indent << "  \"assembly_cache_misses\": "
           << counters.assemblyCacheMisses << ",\n"
           << indent << "  \"assembly_cache_pruned_hits\": "
           << counters.assemblyCachePrunedHits << ",\n"
           << indent << "  \"assembly_cache_updated_hits\": "
           << counters.assemblyCacheUpdatedHits << ",\n"
           << indent << "  \"pair_bound_cache_lookups\": "
           << counters.pairBoundCacheLookups << ",\n"
           << indent << "  \"pair_bound_cache_hits\": "
           << counters.pairBoundCacheHits << ",\n"
           << indent << "  \"pair_bound_cache_misses\": "
           << counters.pairBoundCacheMisses << '\n'
           << indent << '}';
}

inline void writeParallelSearchTelemetry(std::ostream &output)
{
    const ParallelSearchTelemetrySummary &parallel = parallelSearchTelemetry;
    bool uniformLocalThreads = !parallel.localThreadsPerRank.empty();
    uint64_t localThreads = uniformLocalThreads
        ? parallel.localThreadsPerRank.front()
        : 0;
    for (uint64_t count : parallel.localThreadsPerRank)
        if (count != localThreads) uniformLocalThreads = false;

    output << "  \"parallel\": {\n"
           << "    \"enabled\": true,\n"
           << "    \"mode\": \"" << parallel.mode << "\",\n"
           << "    \"aggregation_scope\": \""
           << parallel.aggregationScope << "\",\n"
           << "    \"rank_count\": " << parallel.rankCount << ",\n"
           << "    \"local_threads\": ";
    if (uniformLocalThreads) output << localThreads;
    else output << "null";
    output << ",\n    \"local_threads_per_rank\": [";
    for (size_t i = 0; i < parallel.localThreadsPerRank.size(); ++i)
    {
        if (i != 0) output << ", ";
        output << parallel.localThreadsPerRank[i];
    }
    output << "],\n"
           << "    \"worker_count\": " << parallel.workerCount << ",\n"
           << "    \"elapsed_timing_method\": "
              "\"parallel_region_steady_clock\",\n"
           << "    \"busy_timing_method\": "
              "\"elapsed_minus_scheduler_idle_time\",\n"
           << "    \"branch_scheduler\": {\n"
           << "      \"strategy\": "
              "\"dynamic_leases_with_static_mpi_rank_partition\",\n"
           << "      \"lease_size\": "
           << parallel.branchLeaseSize << ",\n"
           << "      \"rank_partition_count\": "
           << parallel.rankCount << ",\n"
           << "      \"adaptive_splitting\": {\n"
           << "        \"minimum_queued_tasks_per_worker\": "
           << parallelMinimumQueuedTasksPerWorker << ",\n"
           << "        \"target_queued_tasks_per_worker\": "
           << parallelTargetQueuedTasksPerWorker << ",\n"
           << "        \"maximum_queued_tasks_per_worker\": "
           << parallelMaximumQueuedTasksPerWorker << ",\n"
           << "        \"maximum_depth\": "
           << parallelMaximumTaskDepth << ",\n"
           << "        \"warm_start\": \"largest_duplicate_first\"\n"
           << "      },\n"
           << "      \"complete\": "
           << (parallel.branchSchedulerComplete ? "true" : "false") << "\n"
           << "    },\n"
           << "    \"branch_scan_complete\": "
           << (parallel.branchScanComplete ? "true" : "false") << ",\n"
           << "    \"aggregate\": {\n"
           << "      \"branch_candidates\": ";
    if (parallel.branchCandidateCountsAgree)
        output << parallel.branchCandidateCount;
    else output << "null";
    output << ",\n"
           << "      \"branch_leases\": "
           << parallel.branchLeaseCount << ",\n"
           << "      \"branch_assignments\": "
           << parallel.branchAssignmentCount << ",\n"
           << "      \"depth_two_tasks_spawned\": "
           << parallel.depthTwoTaskSpawnCount << ",\n"
           << "      \"depth_two_tasks_executed\": "
           << parallel.depthTwoTaskExecutionCount << ",\n"
           << "      \"deeper_tasks_spawned\": "
           << parallel.deeperTaskSpawnCount << ",\n"
           << "      \"deeper_tasks_executed\": "
           << parallel.deeperTaskExecutionCount << ",\n"
           << "      \"task_steal_attempts\": "
           << parallel.taskStealAttemptCount << ",\n"
           << "      \"task_steals\": "
           << parallel.taskStealCount << ",\n"
           << "      \"local_task_executions\": "
           << parallel.localTaskExecutionCount << ",\n"
           << "      \"scheduler_idle_waits\": "
           << parallel.schedulerIdleWaitCount << ",\n"
           << "      \"scheduler_idle_nanoseconds\": "
           << parallel.schedulerIdleNanoseconds << ",\n"
           << "      \"deep_refill_activations\": "
           << parallel.deepRefillActivationCount << ",\n"
           << "      \"task_queue_high_watermark\": "
           << parallel.taskQueueHighWatermark << ",\n"
           << "      \"maximum_task_depth_executed\": "
           << parallel.maximumTaskDepthExecuted << ",\n"
           << "      \"proactive_tail_refills\": "
           << parallel.proactiveTailRefillCount << ",\n"
           << "      \"warm_start_branches\": "
           << parallel.warmStartBranchCount << ",\n"
           << "      \"elapsed_nanoseconds\": "
           << parallel.elapsedNanoseconds << ",\n"
           << "      \"worker_elapsed_nanoseconds\": "
           << parallel.workerElapsedNanoseconds << ",\n"
           << "      \"worker_busy_nanoseconds\": "
           << parallel.workerBusyNanoseconds << ",\n"
           << "      \"counters\": ";
    writeAllSearchTelemetryCounters(
        output,
        parallel.aggregateCounters,
        "      "
    );
    output << "\n    },\n"
           << "    \"workers\": [\n";
    for (size_t workerIndex = 0; workerIndex < parallel.workers.size();
         ++workerIndex)
    {
        const ParallelSearchWorkerTelemetry &worker =
            parallel.workers[workerIndex];
        output << "      {\n"
               << "        \"rank\": " << worker.mpiRank << ",\n"
               << "        \"local_worker_index\": "
               << worker.localWorkerIndex << ",\n"
               << "        \"global_worker_index\": "
               << worker.globalWorkerIndex << ",\n"
               << "        \"rank_partition\": {\n"
               << "          \"index\": "
               << worker.rankPartitionIndex << ",\n"
               << "          \"count\": "
               << worker.rankPartitionCount << "\n"
               << "        },\n"
               << "        \"branch_candidates\": "
               << worker.branchCandidates << ",\n"
               << "        \"branch_leases\": "
               << worker.branchLeases << ",\n"
               << "        \"branch_assignments\": "
               << worker.branchAssignments << ",\n"
               << "        \"depth_two_tasks_spawned\": "
               << worker.depthTwoTasksSpawned << ",\n"
               << "        \"depth_two_tasks_executed\": "
               << worker.depthTwoTasksExecuted << ",\n"
               << "        \"deeper_tasks_spawned\": "
               << worker.deeperTasksSpawned << ",\n"
               << "        \"deeper_tasks_executed\": "
               << worker.deeperTasksExecuted << ",\n"
               << "        \"task_steal_attempts\": "
               << worker.taskStealAttempts << ",\n"
               << "        \"task_steals\": "
               << worker.taskSteals << ",\n"
               << "        \"local_task_executions\": "
               << worker.localTaskExecutions << ",\n"
               << "        \"scheduler_idle_waits\": "
               << worker.schedulerIdleWaits << ",\n"
               << "        \"scheduler_idle_nanoseconds\": "
               << worker.schedulerIdleNanoseconds << ",\n"
               << "        \"deep_refill_activations\": "
               << worker.deepRefillActivations << ",\n"
               << "        \"task_queue_high_watermark\": "
               << worker.taskQueueHighWatermark << ",\n"
               << "        \"maximum_task_depth_executed\": "
               << worker.maximumTaskDepthExecuted << ",\n"
               << "        \"proactive_tail_refills\": "
               << worker.proactiveTailRefills << ",\n"
               << "        \"warm_start_branches\": "
               << worker.warmStartBranches << ",\n"
               << "        \"elapsed_nanoseconds\": "
               << worker.elapsedNanoseconds << ",\n"
               << "        \"busy_nanoseconds\": "
               << worker.busyNanoseconds << ",\n"
               << "        \"processed_graph\": {\n"
               << "          \"atoms\": " << worker.processedAtoms << ",\n"
               << "          \"edges\": " << worker.processedEdges << ",\n"
               << "          \"active_mask_words\": "
               << worker.activeMaskWords << ",\n"
               << "          \"residual_cache_eligible\": "
               << (worker.residualCacheEligible != 0 ? "true" : "false")
               << "\n        },\n"
               << "        \"phases\": {\n";
        for (size_t i = 0;
             i < static_cast<size_t>(SearchTelemetryPhase::count);
             ++i)
        {
            const SearchTelemetryPhase phase =
                static_cast<SearchTelemetryPhase>(i);
            output << "          \"" << searchTelemetryPhaseName(phase)
                   << "\": {\n"
                   << "            \"wall_nanoseconds\": "
                   << worker.phaseWallNanoseconds[i] << ",\n"
                   << "            \"activations\": "
                   << worker.phaseActivations[i] << "\n"
                   << "          }";
            if (i + 1 < static_cast<size_t>(SearchTelemetryPhase::count))
                output << ',';
            output << '\n';
        }
        output << "        },\n"
               << "        \"counters\": ";
        writeAllSearchTelemetryCounters(output, worker.counters, "        ");
        output << "\n      }";
        if (workerIndex + 1 < parallel.workers.size()) output << ',';
        output << '\n';
    }
    output << "    ]\n"
           << "  }";
}

inline void writeJsonRatio(
    std::ostream &output,
    uint64_t numerator,
    uint64_t denominator
)
{
    if (denominator == 0)
    {
        output << "null";
        return;
    }
    output << static_cast<double>(numerator) /
        static_cast<double>(denominator);
}

inline void writeJsonMemoryValue(
    std::ostream &output,
    bool available,
    uint64_t value
)
{
    if (available) output << value;
    else output << "null";
}

inline bool writeSearchTelemetry(const std::string &filename)
{
    if (!searchTelemetryCompiled || !searchTelemetryEnabled) return true;
    finaliseSearchTelemetry();

    std::ofstream output(filename);
    if (!output.is_open())
    {
        std::cerr << "error: could not open output file '"
                  << filename << "'\n";
        return false;
    }
    const SearchTelemetryCounters &counters = searchTelemetry.counters;
    const uint64_t canonicalLookups =
        counters.canonicalisationMaskCacheHits +
        counters.canonicalisationMaskCacheMisses;
    const uint64_t canonicalClassLookups =
        counters.canonicalClassInsertions + counters.canonicalClassReuses;
    const uint64_t residualLookupCount =
        counters.residualCacheHits + counters.residualCacheMisses;
    const uint64_t assemblyLookupCount =
        counters.assemblyCacheHits + counters.assemblyCacheMisses;
    const uint64_t pairBoundLookupCount =
        counters.pairBoundCacheHits + counters.pairBoundCacheMisses;

    bool anyActivatedPhase = false;
    bool allActivatedPhasePeaksExact = true;
    uint64_t overallResidentPeak = 0;
    for (const SearchTelemetryPhaseStats &phase : searchTelemetry.phases)
    {
        if (phase.activations == 0) continue;
        anyActivatedPhase = true;
        if (!phase.exactResidentPeak)
        {
            allActivatedPhasePeaksExact = false;
            continue;
        }
        overallResidentPeak = std::max(
            overallResidentPeak,
            phase.peakResidentKiB
        );
    }
    const bool overallResidentPeakAvailable =
        anyActivatedPhase && allActivatedPhasePeaksExact;

    output << "{\n"
           << "  \"schema_version\": 1,\n"
           << "  \"processed_graph\": {\n"
           << "    \"atoms\": " << searchTelemetry.processedAtoms << ",\n"
           << "    \"edges\": " << searchTelemetry.processedEdges << ",\n"
           << "    \"active_mask_words\": "
           << searchTelemetry.activeMaskWords << "\n"
           << "  },\n"
           << "  \"counters\": {\n"
           << "    \"retained_mask_attempts\": "
           << counters.retainedMaskAttempts << ",\n"
           << "    \"retained_masks\": " << counters.retainedMasks << ",\n"
           << "    \"duplicate_mask_attempts\": "
           << counters.duplicateMaskAttempts << ",\n"
           << "    \"rejected_masks\": " << counters.rejectedMasks << ",\n"
           << "    \"matching_visits\": " << counters.matchingVisits << ",\n"
           << "    \"canonicalisation_calls\": "
           << counters.canonicalisationCalls << ",\n"
           << "    \"vf2_calls\": " << counters.vf2Calls << ",\n"
           << "    \"vf2_matches\": " << counters.vf2Matches << "\n"
           << "  },\n"
           << "  \"caches\": {\n"
           << "    \"canonical_mask\": {\n"
           << "      \"hits\": "
           << counters.canonicalisationMaskCacheHits << ",\n"
           << "      \"misses\": "
           << counters.canonicalisationMaskCacheMisses << ",\n"
           << "      \"hit_rate\": ";
    writeJsonRatio(
        output,
        counters.canonicalisationMaskCacheHits,
        canonicalLookups
    );
    output << "\n    },\n"
           << "    \"canonical_class\": {\n"
           << "      \"insertions\": "
           << counters.canonicalClassInsertions << ",\n"
           << "      \"reuses\": " << counters.canonicalClassReuses << ",\n"
           << "      \"reuse_rate\": ";
    writeJsonRatio(output, counters.canonicalClassReuses, canonicalClassLookups);
    output << "\n    },\n"
           << "    \"residual_decomposition\": {\n"
           << "      \"eligible_for_processed_graph\": "
           << (searchTelemetry.residualCacheEligible ? "true" : "false")
           << ",\n"
           << "      \"requests\": "
           << counters.residualDecompositionRequests << ",\n"
           << "      \"eligible_requests\": "
           << counters.residualCacheEligibleRequests << ",\n"
           << "      \"small_molecule_bypasses\": "
           << counters.residualCacheSmallMoleculeBypasses << ",\n"
           << "      \"wide_molecule_bypasses\": "
           << counters.residualCacheWideMoleculeBypasses << ",\n"
           << "      \"small_residual_bypasses\": "
           << counters.residualCacheSmallResidualBypasses << ",\n"
           << "      \"first_occurrence_bypasses\": "
           << counters.residualCacheFirstOccurrenceBypasses << ",\n"
           << "      \"runtime_disabled_bypasses\": "
           << counters.residualCacheRuntimeDisabledBypasses << ",\n"
           << "      \"lookups\": " << counters.residualCacheLookups << ",\n"
           << "      \"hits\": " << counters.residualCacheHits << ",\n"
           << "      \"misses\": " << counters.residualCacheMisses << ",\n"
           << "      \"admissions\": "
           << counters.residualCacheAdmissions << ",\n"
           << "      \"lookup_hit_rate\": ";
    writeJsonRatio(output, counters.residualCacheHits, residualLookupCount);
    output << ",\n      \"request_hit_rate\": ";
    writeJsonRatio(
        output,
        counters.residualCacheHits,
        counters.residualDecompositionRequests
    );
    output << "\n    },\n"
           << "    \"assembly_state\": {\n"
           << "      \"lookups\": " << counters.assemblyCacheLookups << ",\n"
           << "      \"hits\": " << counters.assemblyCacheHits << ",\n"
           << "      \"misses\": " << counters.assemblyCacheMisses << ",\n"
           << "      \"pruned_hits\": "
           << counters.assemblyCachePrunedHits << ",\n"
           << "      \"updated_hits\": "
           << counters.assemblyCacheUpdatedHits << ",\n"
           << "      \"hit_rate\": ";
    writeJsonRatio(output, counters.assemblyCacheHits, assemblyLookupCount);
    output << "\n    },\n"
           << "    \"pair_bound\": {\n"
           << "      \"lookups\": " << counters.pairBoundCacheLookups << ",\n"
           << "      \"hits\": " << counters.pairBoundCacheHits << ",\n"
           << "      \"misses\": " << counters.pairBoundCacheMisses << ",\n"
           << "      \"hit_rate\": ";
    writeJsonRatio(output, counters.pairBoundCacheHits, pairBoundLookupCount);
    output << "\n    }\n"
           << "  },\n"
           << "  \"memory\": {\n"
           << "    \"method\": \"";
    if (!searchTelemetry.collectPhaseMemory)
        output << "disabled_parallel";
    else
    {
#ifdef __linux__
        output << "linux_proc_vmhwm_reset";
#else
        output << "unavailable";
#endif
    }
    output << "\",\n"
           << "    \"phase_peaks_are_absolute_not_additive\": true,\n"
           << "    \"phase_peaks_complete\": "
           << (overallResidentPeakAvailable ? "true" : "false") << ",\n"
           << "    \"overall_peak_resident_kib\": ";
    writeJsonMemoryValue(
        output,
        overallResidentPeakAvailable,
        overallResidentPeak
    );
    output << ",\n    \"process_peak_virtual_kib\": ";
    writeJsonMemoryValue(
        output,
        searchTelemetry.finalMemory.virtualAvailable,
        searchTelemetry.finalMemory.virtualPeakKiB
    );
    output << ",\n    \"phases\": {\n";
    for (size_t i = 0; i < static_cast<size_t>(SearchTelemetryPhase::count); i++)
    {
        const SearchTelemetryPhase phase = static_cast<SearchTelemetryPhase>(i);
        const SearchTelemetryPhaseStats &stats = searchTelemetry.phases[i];
        const bool exact = stats.exactResidentPeak;
        output << "      \"" << searchTelemetryPhaseName(phase) << "\": {\n"
               << "        \"clock_ticks\": " << stats.clockTicks << ",\n"
               << "        \"activations\": " << stats.activations << ",\n"
               << "        \"start_rss_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.startResidentAvailable,
            stats.startResidentKiB
        );
        output << ",\n        \"peak_rss_kib\": ";
        writeJsonMemoryValue(output, exact, stats.peakResidentKiB);
        output << ",\n        \"end_rss_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.endResidentAvailable,
            stats.endResidentKiB
        );
        output << ",\n        \"start_virtual_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.startVirtualAvailable,
            stats.startVirtualKiB
        );
        output << ",\n        \"end_virtual_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.endVirtualAvailable,
            stats.endVirtualKiB
        );
        output << "\n      }";
        if (i + 1 < static_cast<size_t>(SearchTelemetryPhase::count))
            output << ',';
        output << '\n';
    }
    output << "    }\n"
           << "  }";
    if (parallelSearchTelemetry.enabled)
    {
        output << ",\n";
        writeParallelSearchTelemetry(output);
    }
    output << "\n}\n";
    output.close();
    if (!output)
    {
        std::cerr << "error: could not write output file '"
                  << filename << "'\n";
        return false;
    }
    return true;
}

#endif
