#define PARALLELASSEMBLYCPP_NO_MAIN
#include "../v5/main.cpp"

#include <cassert>

SearchContext makeTransferContext(size_t edgeCount)
{
    SearchContext context;
    string atomType = "C";
    context.canonicalSeed.atomInterner.emplace(atomType, 1);
    for (size_t atom = 0; atom <= edgeCount; ++atom)
        context.processedMolecule.addAtom(atomType);
    for (size_t edge = 0; edge < edgeCount; ++edge)
        context.processedMolecule.addBond(
            static_cast<int>(edge),
            static_cast<int>(edge + 1),
            1
        );
    context.universeEdges = context.processedMolecule.writeEdgeList();
    context.bondCount = static_cast<unsigned int>(edgeCount);
    context.rootAssemblyIndex = static_cast<int>(edgeCount) - 1;
    context.startedAt = clock();
    return context;
}

assemblyState makeTransferState(size_t edgeCount, int canonicalId = 42)
{
    EdgeMask mask;
    for (size_t edge = 0; edge < edgeCount; ++edge) mask.set(edge);
    assemblyState state;
    state.appendFragment(mask, static_cast<int>(edgeCount), canonicalId, true);
    state.sumDupBonds = 2;
    return state;
}

ParallelTaskScheduler makeTransferScheduler()
{
    return ParallelTaskScheduler(0, 0, 1, 2, 1, false, false, 2);
}

parallelSearchTaskDescriptor takeTransfer(
    ParallelTaskScheduler &scheduler,
    const assemblyState &state,
    size_t receiver = 0
)
{
    assert(scheduler.tryEnqueueTask(0, state, 0, 2));
    parallelSearchTaskDescriptor task;
    assert(scheduler.nextWork(receiver, task) ==
        ParallelTaskScheduler::WorkAvailability::task);
    return task;
}

void testCanonicalValidityAndWideMasks()
{
    SearchContext context = makeTransferContext(130);
    configureParallelWorker(context, 0);
    {
        ParallelTaskScheduler scheduler = makeTransferScheduler();
        assemblyState source = makeTransferState(130, 7);
        parallelSearchTaskDescriptor task = takeTransfer(scheduler, source, 1);
        assert(task.originWorkerIndex == 0);
        assert(task.fragments.front().canonicalId == 7);
        assert(task.fragmentWords.size() == 3);
        assert(task.fragments.front().connected);
        assert(task.sumDupBonds == source.sumDupBonds);
        assemblyState restored;

        reconstructParallelTask(task, restored, 0);
        assert(restored.fragments.front().canonicalId == 7);
        assert(restored.fragments.front().mask == source.fragments.front().mask);
        assert(restored.fragments.front().edgeCount == 130);
        assert(restored.fragments.front().connected);
        assert(restored.sumDupBonds == source.sumDupBonds);

        reconstructParallelTask(task, restored, 1);
        assert(restored.fragments.front().canonicalId == unknownCanonicalId);

        task.canonicalSeedSize = 8;
        reconstructParallelTask(task, restored, 1);
        assert(restored.fragments.front().canonicalId == 7);
        task.canonicalSeedSize = 7;
        reconstructParallelTask(task, restored, 1);
        assert(restored.fragments.front().canonicalId == unknownCanonicalId);

        sharedCanonicalIdRegistry registry(7);
        sharedCanonicalRegistry = &registry;
        reconstructParallelTask(task, restored, 1);
        assert(restored.fragments.front().canonicalId == unknownCanonicalId);
        task.sharedCanonicalIds = true;
        reconstructParallelTask(task, restored, 1);
        assert(restored.fragments.front().canonicalId == 7);
        sharedCanonicalRegistry = nullptr;
        reconstructParallelTask(task, restored, 1);
        assert(restored.fragments.front().canonicalId == unknownCanonicalId);

        task.fragments.front().canonicalId = unknownCanonicalId;
        task.canonicalSeedSize = 8;
        reconstructParallelTask(task, restored, 0);
        assert(restored.fragments.front().canonicalId == unknownCanonicalId);

        task.fragmentWords.pop_back();
        bool rejectedIncompleteMask = false;
        try
        {
            reconstructParallelTask(task, restored, 0);
        }
        catch (const logic_error &)
        {
            rejectedIncompleteMask = true;
        }
        assert(rejectedIncompleteMask);
        scheduler.recycleTask(task);
        scheduler.completeTask(2);
    }
    clearParallelWorkerMasks();
}

void testReusedBuffersAndBounds()
{
    SearchContext context = makeTransferContext(16);
    configureParallelWorker(context, 0);
    {
        assemblyState source = makeTransferState(16);
        {
            ParallelTaskScheduler scheduler = makeTransferScheduler();
            parallelSearchTaskDescriptor first = takeTransfer(scheduler, source);
            const auto *fragments = first.fragments.data();
            const auto *words = first.fragmentWords.data();
            scheduler.recycleTask(first);
            scheduler.completeTask(2);
            parallelSearchTaskDescriptor second = takeTransfer(scheduler, source, 1);
            assert(second.fragments.data() == fragments);
            assert(second.fragmentWords.data() == words);
            scheduler.recycleTask(second);
            scheduler.completeTask(2);

            // Bursts can return more buffers than the per-owner pool keeps.
            for (size_t index = 0;
                 index < ParallelTaskScheduler::maximumRetainedBuffersPerWorker + 4;
                 ++index)
            {
                parallelSearchTaskDescriptor extra;
                extra.fragments.reserve(1);
                extra.fragmentWords.reserve(1);
                scheduler.recycleTask(extra);
            }
            const size_t createdBefore = searchTaskBuffersCreated;
            const size_t reusedBefore = searchTaskBuffersReused;
            vector<parallelSearchTaskDescriptor> held;
            for (size_t index = 0;
                 index <= ParallelTaskScheduler::maximumRetainedBuffersPerWorker;
                 ++index)
            {
                held.push_back(takeTransfer(scheduler, source));
                scheduler.completeTask(2);
            }
            assert(searchTaskBuffersReused - reusedBefore ==
                ParallelTaskScheduler::maximumRetainedBuffersPerWorker);
            assert(searchTaskBuffersCreated - createdBefore == 1);
        }
        {
            ParallelTaskScheduler scheduler = makeTransferScheduler();
            parallelSearchTaskDescriptor held = takeTransfer(scheduler, source);
            scheduler.completeTask(2);
            parallelSearchTaskDescriptor oversized;
            oversized.fragmentWords.reserve(
                ParallelTaskScheduler::maximumRetainedBufferBytesPerWorker /
                    sizeof(uint64_t) + 1
            );
            scheduler.recycleTask(oversized);
            assert(oversized.fragmentWords.capacity() == 0);

            // Each of these fits alone, but two exceed the byte allowance.
            for (size_t index = 0; index < 3; ++index)
            {
                parallelSearchTaskDescriptor extra;
                extra.fragments.reserve(1);
                extra.fragmentWords.reserve(
                    ParallelTaskScheduler::maximumRetainedBufferBytesPerWorker /
                        (2 * sizeof(uint64_t)) + 1
                );
                scheduler.recycleTask(extra);
            }
            const size_t reusedBefore = searchTaskBuffersReused;
            parallelSearchTaskDescriptor first = takeTransfer(scheduler, source);
            scheduler.completeTask(2);
            parallelSearchTaskDescriptor second = takeTransfer(scheduler, source);
            scheduler.completeTask(2);
            assert(searchTaskBuffersReused - reusedBefore == 1);
            scheduler.cancel();
            scheduler.recycleTask(first);
            scheduler.recycleTask(second);
            assert(!scheduler.tryEnqueueTask(0, source, 0, 2));
        }
    }
    clearParallelWorkerMasks();
}

void testExecutionPruningCancellationAndException()
{
    SearchContext context = makeTransferContext(16);
    // A leaf-only DAG enters the recursive search and then exhausts its work.
    context.dag.resize(1);
    context.dag.front().nodes.resize(context.universeEdges.size());
    configureParallelWorker(context, 0);
    {
        WorkerContext worker(context);
        assemblyState source = makeTransferState(8, unknownCanonicalId);
        source.fragments.front().canonicalId = canonise(source.fragments.front().mask);
        bitsetHashTable.clear();

        ParallelTaskScheduler scheduler = makeTransferScheduler();
        parallelSearchTaskDescriptor task = takeTransfer(scheduler, source);
        bool immediatelyPruned = true;
        assert(runParallelTask<false>(context, task, worker, immediatelyPruned));
        assert(!immediatelyPruned);
        // Preserved IDs avoid recanonicalising known fragments on restoration.
        assert(bitsetHashTable.empty());
        const size_t prunedBefore = searchTasksImmediatelyPruned;
        assert(runAndRecycleParallelTask<false>(context, task, worker, scheduler));
        assert(searchTasksImmediatelyPruned == prunedBefore + 1);
        assert(scheduler.nextWork(0, task) ==
            ParallelTaskScheduler::WorkAvailability::complete);

        task = takeTransfer(scheduler, source);
        task.lowerBoundAssemblyIndex = worker.assemblyIndex;
        assert(runAndRecycleParallelTask<false>(context, task, worker, scheduler));
        assert(searchTasksImmediatelyPruned == prunedBefore + 2);

        task = takeTransfer(scheduler, source);
        searchCancellationFlag.store(true);
        assert(!runAndRecycleParallelTask<false>(context, task, worker, scheduler));
        searchCancellationFlag.store(false);
        assert(searchTasksImmediatelyPruned == prunedBefore + 2);
        assert(scheduler.nextWork(0, task) ==
            ParallelTaskScheduler::WorkAvailability::complete);

        task = takeTransfer(scheduler, source);
        const auto *words = task.fragmentWords.data();
        task.fragmentWords.clear();
        bool threw = false;
        try
        {
            static_cast<void>(runAndRecycleParallelTask<false>(
                context, task, worker, scheduler
            ));
        }
        catch (const logic_error &)
        {
            threw = true;
        }
        assert(threw);
        assert(searchTasksImmediatelyPruned == prunedBefore + 2);
        assert(scheduler.nextWork(0, task) ==
            ParallelTaskScheduler::WorkAvailability::complete);
        task = takeTransfer(scheduler, source);
        assert(task.fragmentWords.data() == words);
        scheduler.recycleTask(task);
        scheduler.completeTask(2);
        assert(searchTaskExecutionNanoseconds > 0);
        assert(searchTaskSerializationNanoseconds > 0);
    }
    clearParallelWorkerMasks();
}

void testMeasuredMinimumTaskSize()
{
    SearchContext context = makeTransferContext(16);
    configureParallelWorker(context, 0);
    {
        ParallelTaskScheduler scheduler = makeTransferScheduler();
        assemblyState small = makeTransferState(8);
        assemblyState large = makeTransferState(16);
        assert(ParallelTaskScheduler::estimateTaskWork(small) == 28);
        assert(ParallelTaskScheduler::estimateTaskWork(large) == 120);
        parallelSearchTaskDescriptor sample;
        sample.estimatedWorkUnits = 100;
        sample.serializationNanoseconds = 100;
        for (size_t index = 0;
             index < ParallelTaskScheduler::taskCalibrationSamples;
             ++index)
            scheduler.recordTaskExecution(sample, 1000, false);

        // Eightfold amortisation makes the measured minimum 81 edge pairs.
        assert(!scheduler.tryEnqueueTask(0, small, 0, 2));
        parallelSearchTaskDescriptor task = takeTransfer(scheduler, large);
        scheduler.recycleTask(task);
        scheduler.completeTask(2);

        // Losing half the donations to immediate pruning doubles the useful
        // work needed to cover the same serialization cost.
        searchTaskMinimumWorkUnits = 0;
        for (size_t index = 0;
             index < ParallelTaskScheduler::taskCalibrationSamples;
             ++index)
            scheduler.recordTaskExecution(
                sample,
                1000,
                index < ParallelTaskScheduler::taskCalibrationSamples / 2
            );
        assert(searchTaskMinimumWorkUnits == 161);
        assert(!scheduler.tryEnqueueTask(0, large, 0, 2));

        for (size_t index = 0;
             index < ParallelTaskScheduler::taskCalibrationSamples;
             ++index)
            scheduler.recordTaskExecution(sample, 1000, true);
        assert(!scheduler.tryEnqueueTask(0, large, 0, 2));
        size_t probes = 0;
        for (size_t index = 0;
             index < ParallelTaskScheduler::taskCalibrationSamples;
             ++index)
        {
            if (scheduler.tryEnqueueTask(0, small, 0, 2))
            {
                ++probes;
                assert(scheduler.nextWork(0, task) ==
                    ParallelTaskScheduler::WorkAvailability::task);
                scheduler.recycleTask(task);
                scheduler.completeTask(2);
            }
        }
        assert(probes == 1);

        // Useful later work lowers the minimum again; pruning cannot disable
        // donation permanently after the incumbent or available work changes.
        for (size_t index = 0;
             index < ParallelTaskScheduler::taskCalibrationSamples;
             ++index)
            scheduler.recordTaskExecution(sample, 1000000, false);
        task = takeTransfer(scheduler, small);
        scheduler.recycleTask(task);
        scheduler.completeTask(2);
    }
    clearParallelWorkerMasks();
}

#ifdef PARALLELASSEMBLYCPP_USE_OPENMP
void testWideMasksCrossWorkerArenas()
{
    SearchContext context = makeTransferContext(130);
    ParallelTaskScheduler scheduler = makeTransferScheduler();
    omp_set_dynamic(0);
    #pragma omp parallel num_threads(2)
    {
        const size_t index = static_cast<size_t>(omp_get_thread_num());
        assert(omp_get_num_threads() == 2);
        configureParallelWorker(context, index);
        if (index == 0)
        {
            assemblyState source = makeTransferState(130, 7);
            assert(scheduler.tryEnqueueTask(0, source, 0, 2));
        }
        #pragma omp barrier
        if (index == 1)
        {
            parallelSearchTaskDescriptor task;
            assert(scheduler.nextWork(1, task) ==
                ParallelTaskScheduler::WorkAvailability::task);
            assemblyState restored;
            reconstructParallelTask(task, restored, 1);
            assert(restored.fragments.front().mask.count() == 130);
            assert(restored.fragments.front().mask.test(129));
            assert(restored.fragments.front().canonicalId == unknownCanonicalId);
            scheduler.recycleTask(task);
            scheduler.completeTask(2);
        }
        clearParallelWorkerMasks();
    }
    // The consumer returned ordinary vectors to the producer's pool, and both
    // worker mask arenas can be reconfigured after their local masks die.
    configureParallelWorker(context, 0);
    {
        assemblyState source = makeTransferState(130, 9);
        const size_t reusedBefore = searchTaskBuffersReused;
        parallelSearchTaskDescriptor task = takeTransfer(scheduler, source);
        assert(searchTaskBuffersReused == reusedBefore + 1);
        scheduler.recycleTask(task);
        scheduler.completeTask(2);
    }
    clearParallelWorkerMasks();
}
#endif

int main()
{
    suppressSearchOutput = true;
    maximumRuntimeTicks = numeric_limits<unsigned long long>::max();
    testCanonicalValidityAndWideMasks();
    testReusedBuffersAndBounds();
    testExecutionPruningCancellationAndException();
    testMeasuredMinimumTaskSize();
#ifdef PARALLELASSEMBLYCPP_USE_OPENMP
    testWideMasksCrossWorkerArenas();
#endif
    return 0;
}
