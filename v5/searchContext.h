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

/** Row-major aggregate masks whose backing allocations survive frame reuse. */
struct flatEdgeMaskAccumulatorTable
{
    EdgeMaskAccumulatorBuffer masks;
    size_t rows = 0;
    size_t columns = 0;

    void reset(size_t columnCount)
    {
        masks.clear();
        rows = 0;
        columns = columnCount;
    }

    span<EdgeMaskAccumulator> appendRow()
    {
        if (
            columns != 0 &&
            rows == numeric_limits<size_t>::max() / columns
        ) throw length_error("aggregate mask table exceeds capacity");
        const size_t offset = rows * columns;
        masks.resize(offset + columns);
        ++rows;
        return masks.span().subspan(offset, columns);
    }

    span<EdgeMaskAccumulator> row(size_t index)
    {
        if (index >= rows) throw out_of_range("aggregate mask row");
        return masks.span().subspan(index * columns, columns);
    }

    span<const EdgeMaskAccumulator> row(size_t index) const
    {
        if (index >= rows) throw out_of_range("aggregate mask row");
        return masks.span().subspan(index * columns, columns);
    }

    size_t maskCount(size_t rowIndex, size_t columnIndex) const
    {
        if (rowIndex >= rows || columnIndex >= columns)
            throw out_of_range("aggregate mask cell");
        return masks[rowIndex * columns + columnIndex].count();
    }
};

/** Mutable storage reused by one worker at one active recursive depth. */
struct dagAssemblySearchFrame
{
    vector<dagDuplicateClassLevel> duplicateLevels;
    size_t duplicateLevelCount = 0;
    flatEdgeMaskAccumulatorTable targetMasks;
    EdgeMaskAccumulatorBuffer aggregateMasks;
    vi fragSizeListMax;
    vi unrestrictedParentTotals;
    vi pairGenericBoundCache;
    assemblyState candidate;

    void reset(size_t fragmentCount)
    {
        for (dagDuplicateClassLevel &level : duplicateLevels)
            level.reset(fragmentCount);
        duplicateLevelCount = 0;
        targetMasks.reset(fragmentCount);
        aggregateMasks.reset(fragmentCount + 2);
        fragSizeListMax.clear();
        unrestrictedParentTotals.clear();
        pairGenericBoundCache.clear();
        candidate.clearFragments();
        candidate.reserveFragments(fragmentCount + 2);
    }

    dagDuplicateClassLevel &appendDuplicateLevel(size_t fragmentCount)
    {
        if (duplicateLevelCount == duplicateLevels.size())
        {
            duplicateLevels.emplace_back();
            duplicateLevels.back().reset(fragmentCount);
        }
        dagDuplicateClassLevel &level =
            duplicateLevels[duplicateLevelCount++];
        return level;
    }
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
    deque<dagAssemblySearchFrame> recursiveFrames;
    size_t recursiveDepth = 0;
    // Parallel workers use these to route donated descendants back to their
    // owner deque and to derive an absolute, bounded search-tree depth. Serial
    // searches retain the harmless defaults.
    size_t parallelWorkerIndex = 0;
    unsigned int parallelTaskDepth = 1;

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

    dagAssemblySearchFrame &acquireRecursiveFrame(size_t fragmentCount)
    {
        if (recursiveDepth == recursiveFrames.size())
            recursiveFrames.emplace_back();
        dagAssemblySearchFrame &frame = recursiveFrames[recursiveDepth++];
        frame.reset(fragmentCount);
        return frame;
    }

    void releaseRecursiveFrame() noexcept
    {
        if (recursiveDepth == 0) terminate();
        --recursiveDepth;
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

/** Balance frame depth on every recursive exit, including cancellation. */
struct dagAssemblySearchFrameScope
{
    assemblySearchStorage &storage;
    dagAssemblySearchFrame &frame;

    dagAssemblySearchFrameScope(
        assemblySearchStorage &_storage,
        size_t fragmentCount
    ):
        storage(_storage),
        frame(storage.acquireRecursiveFrame(fragmentCount)) {}

    ~dagAssemblySearchFrameScope()
    {
        storage.releaseRecursiveFrame();
    }

    dagAssemblySearchFrameScope(const dagAssemblySearchFrameScope &) = delete;
    dagAssemblySearchFrameScope &operator=(
        const dagAssemblySearchFrameScope &
    ) = delete;
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

    explicit WorkerContext(
        const SearchContext &context,
        std::size_t workerIndex = 0
    ):
        fragmentation(
            context.processedMolecule.mg.size(),
            context.universeEdges.size()
        ),
        assemblyIndex(context.rootAssemblyIndex)
    {
        search.parallelWorkerIndex = workerIndex;
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

/** One fragment of a transferable search task, stored without an EdgeMask. */
struct parallelTaskFragmentDescriptor
{
    std::size_t wordOffset = 0;
    std::uint32_t edgeCount = 0;
    bool connected = false;
};

/**
 * An assembly state that may safely move between worker threads.
 *
 * Canonical IDs are deliberately omitted: after the shared root seed they are
 * process-global only when L2 reuse is enabled, and worker-local otherwise.
 * The receiving worker reconstructs the masks and canonicalises them through
 * the active mode before continuing the search.
 */
struct parallelSearchTaskDescriptor
{
    std::vector<parallelTaskFragmentDescriptor> fragments;
    std::vector<std::uint64_t> fragmentWords;
    int sumDupBonds = 0;
    int lowerBoundAssemblyIndex = 0;
    unsigned int depth = 2;
};

inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchDepthTwoTasksSpawned = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchDepthTwoTasksExecuted = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchDeeperTasksSpawned = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchDeeperTasksExecuted = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchTaskStealAttempts = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchTaskSteals = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchLocalTaskExecutions = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchSchedulerIdleWaits = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::uint64_t searchSchedulerIdleNanoseconds = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchDeepRefillActivations = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchTaskQueueHighWatermark = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL unsigned int searchMaximumTaskDepthExecuted = 0;
inline ASSEMBLYCPP_SEARCH_LOCAL std::size_t searchWarmStartBranches = 0;

/** Keep independent scheduler ownership domains off the same cache line. */
template<typename Value>
struct alignas(schedulerCacheLineBytes) schedulerAtomicValue
{
    std::atomic<Value> value;

    explicit schedulerAtomicValue(Value initial = Value{}) noexcept:
        value(initial) {}

    schedulerAtomicValue(const schedulerAtomicValue &) = delete;
    schedulerAtomicValue &operator=(const schedulerAtomicValue &) = delete;
};

static_assert(
    alignof(schedulerAtomicValue<std::size_t>) >= schedulerCacheLineBytes
);
static_assert(
    sizeof(schedulerAtomicValue<std::size_t>) >= schedulerCacheLineBytes
);

struct alignas(schedulerCacheLineBytes) schedulerTaskDepthCounts
{
    std::array<
        std::atomic<std::size_t>,
        parallelMaximumTaskDepth + 1
    > values;
};

static_assert(sizeof(schedulerTaskDepthCounts) >= schedulerCacheLineBytes);

struct alignas(schedulerCacheLineBytes) parallelWorkerTaskDeque
{
    std::atomic<std::size_t> ready{0};
    std::mutex mutex;
    // Keep regular searches allocation-free; construct a deque only after
    // live starvation or the existing shallow policy requests donation.
    std::optional<std::deque<parallelSearchTaskDescriptor>> tasks;
    // Only this deque's owner updates its rotating steal origin.
    std::size_t stealOffset = 0;
};

/**
 * Rank-local adaptive work controller.
 *
 * Root jobs retain their stable modulo MPI partition and are claimed in
 * leases. Shallow donation retains the existing sparse/proactive policy.
 * Deeper donation is armed one level at a time only after idle workers have
 * failed to find local or stealable work. Each worker owns a LIFO deque; other
 * workers steal older tasks from its opposite end. Queue size and task depth
 * are both bounded.
 */
class ParallelTaskScheduler
{
public:
    enum class WorkAvailability
    {
        task,
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
        taskDeques(std::make_unique<parallelWorkerTaskDeque[]>(
            std::max<std::size_t>(1, localWorkerCount)
        )),
        rootLeaseCursor(0),
        rootJobsOutstanding(rankRootJobCount),
        taskJobsOutstanding(0),
        readyTaskCount(0),
        taskSlotsInUse(0),
        idleWorkers(0),
        requestedDonationDepth(0),
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
        taskTransfersEnabled(localWorkerCount > 1),
        proactiveTailRefillEnabled(
            taskTransfersEnabled &&
            rankRootJobCount >= saturatedProduct(
                targetWatermark,
                parallelPromisingFrontierLeaseSize
            )
        ),
        idleWorkerTrigger(
            localWorkerCount <= 1
                ? 1
                : localWorkerCount < 8
                    ? (localWorkerCount + 1) / 2
                    : std::max<std::size_t>(
                        2,
                        (localWorkerCount + 3) / 4
                    )
        )
    {
        for (auto &outstanding : taskDepthOutstanding.values)
            outstanding.store(0, std::memory_order_relaxed);
        const std::size_t sparseThreshold = std::max<std::size_t>(
            1,
            lowWatermark / 2
        );
        if (taskTransfersEnabled && rankRootJobCount < sparseThreshold)
            requestedDonationDepth.value.store(2, std::memory_order_relaxed);
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
        begin = rootLeaseCursor.value.load(std::memory_order_relaxed);
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
        while (!rootLeaseCursor.value.compare_exchange_weak(
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
            if (!cancelled.load(std::memory_order_relaxed))
            {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                ++searchProactiveTailRefills;
#endif
                requestDonationDepth(2);
                notifyWaiters(true);
            }
        }
        return true;
    }

    [[nodiscard]] bool taskDonationRequested(
        unsigned int taskDepth
    ) const noexcept
    {
        return taskDepth <= requestedDonationDepth.value.load(
            std::memory_order_relaxed
        );
    }

    [[nodiscard]] bool transferableTasksEnabled() const noexcept
    {
        return taskTransfersEnabled;
    }

    /**
     * Queue a bounded descendant on the producer's own deque.
     */
    bool tryEnqueueTask(
        std::size_t workerIndex,
        const assemblyState &state,
        int lowerBoundAssemblyIndex,
        unsigned int taskDepth
    )
    {
        if (
            workerIndex >= workerCount ||
            taskDepth < 2 ||
            taskDepth > maximumTaskDepth ||
            !taskDonationRequested(taskDepth) ||
            cancelled.load(std::memory_order_relaxed) ||
            !reserveTaskSlot(taskDepth)
        ) return false;

        parallelSearchTaskDescriptor task;
        try
        {
            task.sumDupBonds = state.sumDupBonds;
            task.lowerBoundAssemblyIndex = lowerBoundAssemblyIndex;
            task.depth = taskDepth;
            task.fragments.reserve(state.fragments.size());
            const std::size_t wordCount = EdgeMask::activeWordCount();
            if (
                wordCount != 0 &&
                state.fragments.size() >
                    std::numeric_limits<std::size_t>::max() / wordCount
            ) throw std::length_error("parallel task masks exceed capacity");
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
            releaseReservedTask(taskDepth);
            notifyWaiters(true);
            throw;
        }

        if (cancelled.load(std::memory_order_relaxed))
        {
            releaseReservedTask(taskDepth);
            notifyWaiters(true);
            return false;
        }

        parallelWorkerTaskDeque &owner = taskDeques[workerIndex];
        std::size_t ownerQueueSize = 0;
        try
        {
            std::lock_guard<std::mutex> lock(owner.mutex);
            if (!owner.tasks.has_value()) owner.tasks.emplace();
            owner.tasks->push_back(std::move(task));
            ownerQueueSize = owner.tasks->size();
            owner.ready.fetch_add(1, std::memory_order_release);
            readyTaskCount.value.fetch_add(1, std::memory_order_release);
        }
        catch (...)
        {
            releaseReservedTask(taskDepth);
            notifyWaiters(true);
            throw;
        }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        searchTaskQueueHighWatermark = std::max(
            searchTaskQueueHighWatermark,
            ownerQueueSize
        );
#endif
        disableDonationAtTarget();
        notifyWaiters(false);
        return true;
    }

    WorkAvailability nextWork(
        std::size_t workerIndex,
        parallelSearchTaskDescriptor &task
    )
    {
        if (workerIndex >= workerCount)
            throw std::out_of_range("parallel worker deque index");

        if (popLocalTask(workerIndex, task))
        {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            ++searchLocalTaskExecutions;
#endif
            return WorkAvailability::task;
        }

        const std::size_t victimCount = workerCount - 1;
        const std::size_t firstRelativeVictim = victimCount == 0
            ? 0
            : 1 + taskDeques[workerIndex].stealOffset % victimCount;
        for (std::size_t offset = 0; offset < victimCount; ++offset)
        {
            if (readyTaskCount.value.load(std::memory_order_acquire) == 0)
                break;
            const std::size_t relativeVictim = 1 + (
                firstRelativeVictim - 1 + offset
            ) % victimCount;
            const std::size_t victim = (workerIndex + relativeVictim) %
                workerCount;
            if (taskDeques[victim].ready.load(
                std::memory_order_acquire
            ) == 0) continue;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            ++searchTaskStealAttempts;
#endif
            if (stealTask(victim, task))
            {
                taskDeques[workerIndex].stealOffset = relativeVictim;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                ++searchTaskSteals;
#endif
                return WorkAvailability::task;
            }
        }
        if (victimCount != 0)
        {
            taskDeques[workerIndex].stealOffset =
                (taskDeques[workerIndex].stealOffset + 1) % victimCount;
        }
        if (
            rootJobsOutstanding.value.load(std::memory_order_acquire) == 0 &&
            taskJobsOutstanding.value.load(std::memory_order_acquire) == 0
        ) return WorkAvailability::complete;
        return WorkAvailability::wait;
    }

    void completeRootLease(std::size_t completedRootJobs)
    {
        if (!taskTransfersEnabled) return;
        if (completedRootJobs != 0)
        {
            const std::size_t previous = rootJobsOutstanding.value.fetch_sub(
                completedRootJobs,
                std::memory_order_acq_rel
            );
            if (previous < completedRootJobs)
            {
                rootJobsOutstanding.value.fetch_add(
                    completedRootJobs,
                    std::memory_order_relaxed
                );
                throw std::logic_error(
                    "adaptive scheduler completed extra root jobs"
                );
            }
            if (previous == completedRootJobs) notifyWaiters(true);
        }
    }

    void completeTask(unsigned int taskDepth)
    {
        if (taskDepth < 2 || taskDepth > maximumTaskDepth)
            throw std::logic_error("adaptive scheduler task depth is invalid");
        const std::size_t previousDepth =
            taskDepthOutstanding.values[taskDepth]
            .fetch_sub(1, std::memory_order_acq_rel);
        if (previousDepth == 0)
        {
            taskDepthOutstanding.values[taskDepth].fetch_add(
                1,
                std::memory_order_relaxed
            );
            throw std::logic_error(
                "adaptive scheduler completed an extra task at this depth"
            );
        }
        const std::size_t previous = taskJobsOutstanding.value.fetch_sub(
            1,
            std::memory_order_acq_rel
        );
        if (previous == 0)
        {
            taskJobsOutstanding.value.fetch_add(1, std::memory_order_relaxed);
            taskDepthOutstanding.values[taskDepth].fetch_add(
                1,
                std::memory_order_relaxed
            );
            throw std::logic_error("adaptive scheduler completed an extra task");
        }
        if (previous == 1) notifyWaiters(true);
    }

    void waitForWork(std::size_t workerIndex)
    {
        if (workerIndex >= workerCount)
            throw std::out_of_range("parallel worker wait index");
        std::unique_lock<std::mutex> lock(waitMutex);
        idleWorkers.value.fetch_add(1, std::memory_order_acq_rel);
        requestStarvationRefill();
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        ++searchSchedulerIdleWaits;
        const auto idleStarted = std::chrono::steady_clock::now();
#endif
        // The timeout only polls process signals. Refill has already been
        // requested from live starvation telemetry before sleeping.
        condition.wait_for(
            lock,
            std::chrono::milliseconds(1),
            [&]
            {
                return
                    cancelled.load(std::memory_order_relaxed) ||
                    readyTaskCount.value.load(std::memory_order_acquire) != 0 ||
                    rootLeaseCursor.value.load(std::memory_order_relaxed) <
                        rankRootJobCount ||
                    (rootJobsOutstanding.value.load(
                        std::memory_order_acquire
                    ) == 0 && taskJobsOutstanding.value.load(
                        std::memory_order_acquire
                    ) == 0);
            }
        );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        const auto idleElapsed = std::chrono::duration_cast<
            std::chrono::nanoseconds
        >(std::chrono::steady_clock::now() - idleStarted).count();
        if (idleElapsed > 0)
            searchSchedulerIdleNanoseconds += static_cast<std::uint64_t>(
                idleElapsed
            );
#endif
        idleWorkers.value.fetch_sub(1, std::memory_order_acq_rel);
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

    [[nodiscard]] std::size_t queuedTaskEstimate() const noexcept
    {
        const std::size_t claimed = std::min(
            rootLeaseCursor.value.load(std::memory_order_relaxed),
            rankRootJobCount
        );
        const std::size_t remainingRootJobs = rankRootJobCount - claimed;
        const std::size_t outstandingRoots = rootJobsOutstanding.value.load(
            std::memory_order_relaxed
        );
        const std::size_t claimedOutstanding =
            outstandingRoots > remainingRootJobs
            ? outstandingRoots - remainingRootJobs
            : 0;
        const std::size_t activeRoots = std::min(
            claimedOutstanding,
            workerCount
        );
        const std::size_t slots = taskSlotsInUse.value.load(
            std::memory_order_relaxed
        );
        const std::size_t outstandingTasks = taskJobsOutstanding.value.load(
            std::memory_order_relaxed
        );
        const std::size_t activeTasks = outstandingTasks > slots
            ? std::min(outstandingTasks - slots, workerCount)
            : 0;
        const std::size_t activeEstimate = std::min(
            workerCount,
            activeRoots > workerCount - std::min(activeTasks, workerCount)
                ? workerCount
                : activeRoots + activeTasks
        );
        if (activeEstimate >
            std::numeric_limits<std::size_t>::max() - remainingRootJobs)
        {
            return std::numeric_limits<std::size_t>::max();
        }
        const std::size_t roots = remainingRootJobs + activeEstimate;
        if (slots > std::numeric_limits<std::size_t>::max() - roots)
            return std::numeric_limits<std::size_t>::max();
        return roots + slots;
    }

    bool reserveTaskSlot(unsigned int taskDepth) noexcept
    {
        std::size_t slots = taskSlotsInUse.value.load(
            std::memory_order_relaxed
        );
        while (true)
        {
            if (
                slots >= maximumQueuedTasks ||
                queuedTaskEstimate() >= targetWatermark ||
                !taskDonationRequested(taskDepth) ||
                cancelled.load(std::memory_order_relaxed)
            ) return false;
            if (taskSlotsInUse.value.compare_exchange_weak(
                slots,
                slots + 1,
                std::memory_order_acq_rel,
                std::memory_order_relaxed
            )) break;
        }
        taskJobsOutstanding.value.fetch_add(1, std::memory_order_acq_rel);
        taskDepthOutstanding.values[taskDepth].fetch_add(
            1,
            std::memory_order_acq_rel
        );
        return true;
    }

    void releaseReservedTask(unsigned int taskDepth) noexcept
    {
        taskDepthOutstanding.values[taskDepth].fetch_sub(
            1,
            std::memory_order_acq_rel
        );
        taskJobsOutstanding.value.fetch_sub(1, std::memory_order_acq_rel);
        taskSlotsInUse.value.fetch_sub(1, std::memory_order_acq_rel);
    }

    void consumePublishedTaskSlot() noexcept
    {
        const std::size_t ready = readyTaskCount.value.fetch_sub(
            1,
            std::memory_order_acq_rel
        );
        const std::size_t slots = taskSlotsInUse.value.fetch_sub(
            1,
            std::memory_order_acq_rel
        );
        if (ready == 0 || slots == 0) std::terminate();
    }

    bool popLocalTask(
        std::size_t workerIndex,
        parallelSearchTaskDescriptor &task
    )
    {
        if (readyTaskCount.value.load(std::memory_order_acquire) == 0)
            return false;
        parallelWorkerTaskDeque &owner = taskDeques[workerIndex];
        std::lock_guard<std::mutex> lock(owner.mutex);
        if (!owner.tasks.has_value() || owner.tasks->empty()) return false;
        task = std::move(owner.tasks->back());
        owner.tasks->pop_back();
        if (owner.ready.fetch_sub(1, std::memory_order_acq_rel) == 0)
            std::terminate();
        consumePublishedTaskSlot();
        return true;
    }

    bool stealTask(
        std::size_t victimIndex,
        parallelSearchTaskDescriptor &task
    )
    {
        if (
            readyTaskCount.value.load(std::memory_order_acquire) == 0 ||
            taskDeques[victimIndex].ready.load(std::memory_order_acquire) == 0
        )
            return false;
        parallelWorkerTaskDeque &victim = taskDeques[victimIndex];
        std::unique_lock<std::mutex> lock(victim.mutex, std::try_to_lock);
        if (
            !lock.owns_lock() ||
            !victim.tasks.has_value() ||
            victim.tasks->empty()
        ) return false;
        task = std::move(victim.tasks->front());
        victim.tasks->pop_front();
        if (victim.ready.fetch_sub(1, std::memory_order_acq_rel) == 0)
            std::terminate();
        consumePublishedTaskSlot();
        return true;
    }

    bool requestDonationDepth(unsigned int depth) noexcept
    {
        depth = std::min(depth, maximumTaskDepth);
        unsigned int requested = requestedDonationDepth.value.load(
            std::memory_order_relaxed
        );
        while (requested < depth)
        {
            if (requestedDonationDepth.value.compare_exchange_weak(
                requested,
                depth,
                std::memory_order_release,
                std::memory_order_relaxed
            )) return true;
        }
        return false;
    }

    void disableDonationAtTarget() noexcept
    {
        if (queuedTaskEstimate() < targetWatermark) return;
        unsigned int requested = requestedDonationDepth.value.load(
            std::memory_order_relaxed
        );
        while (
            requested != 0 &&
            !requestedDonationDepth.value.compare_exchange_weak(
                requested,
                0,
                std::memory_order_release,
                std::memory_order_relaxed
            )
        ) {}
    }

    void requestStarvationRefill() noexcept
    {
        if (
            !taskTransfersEnabled ||
            cancelled.load(std::memory_order_relaxed) ||
            readyTaskCount.value.load(std::memory_order_acquire) != 0 ||
            idleWorkers.value.load(std::memory_order_acquire) <
                idleWorkerTrigger ||
            queuedTaskEstimate() >= lowWatermark
        ) return;

        unsigned int desiredDepth = 0;
        for (unsigned int depth = maximumTaskDepth; depth > 2; --depth)
        {
            if (taskDepthOutstanding.values[depth - 1].load(
                std::memory_order_acquire
            ) != 0)
            {
                desiredDepth = depth;
                break;
            }
        }
        if (
            desiredDepth == 0 &&
            rootJobsOutstanding.value.load(std::memory_order_acquire) != 0
        ) desiredDepth = 2;
        if (desiredDepth == 0) return;

#ifdef ASSEMBLY_ENABLE_TELEMETRY
        const bool activated = requestDonationDepth(desiredDepth);
        if (desiredDepth > 2 && activated)
            ++searchDeepRefillActivations;
#else
        static_cast<void>(requestDonationDepth(desiredDepth));
#endif
    }

    void notifyWaiters(bool all)
    {
        // Publication happens before acquiring this mutex. A waiter cannot
        // miss the subsequent notification between its predicate and sleep.
        {
            std::lock_guard<std::mutex> lock(waitMutex);
        }
        if (all) condition.notify_all();
        else condition.notify_one();
    }

    const std::size_t globalRootJobCount;
    const std::size_t rankIndex;
    const std::size_t rankCount;
    const std::size_t leaseSize;
    const std::size_t rankRootJobCount;
    const std::size_t workerCount;
    std::unique_ptr<parallelWorkerTaskDeque[]> taskDeques;
    schedulerAtomicValue<std::size_t> rootLeaseCursor;
    schedulerAtomicValue<std::size_t> rootJobsOutstanding;
    schedulerAtomicValue<std::size_t> taskJobsOutstanding;
    schedulerAtomicValue<std::size_t> readyTaskCount;
    schedulerAtomicValue<std::size_t> taskSlotsInUse;
    schedulerAtomicValue<std::size_t> idleWorkers;
    schedulerAtomicValue<unsigned int> requestedDonationDepth;
    schedulerTaskDepthCounts taskDepthOutstanding;
    mutable std::mutex waitMutex;
    std::condition_variable condition;
    const std::size_t lowWatermark;
    const std::size_t targetWatermark;
    const std::size_t maximumQueuedTasks;
    const bool adaptiveRootLeases;
    const bool taskTransfersEnabled;
    const bool proactiveTailRefillEnabled;
    const std::size_t idleWorkerTrigger;
    std::atomic_bool rootFrontierRefillActivated{false};
    std::atomic_bool cancelled{false};
};

ASSEMBLYCPP_SEARCH_LOCAL ParallelTaskScheduler *parallelTaskScheduler = nullptr;
