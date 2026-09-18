#define PARALLELASSEMBLYCPP_NO_MAIN
#include "../src/main.cpp"

#include <cassert>

int main()
{
    suppressSearchOutput = true;
    maximumRuntimeTicks = numeric_limits<unsigned long long>::max();
    molGraph molecule;
    string carbon = "C";
    for (int atom = 0; atom < 33; ++atom) molecule.addAtom(carbon);
    for (int atom = 0; atom < 32; ++atom) molecule.addBond(atom, atom + 1, 1);
    SearchContext context;
    prepareParallelSearchContext(molecule, context);
    assert(context.rootJobs.size() > 1);
    configureParallelWorker(context, 0);
    // Explicit nested pairs exercise a multi-step feasible path independently
    // of homogeneous-path equivalence, which may keep different root pairs.
    context.rootJobs.clear();
    context.rootOccurrences.clear();
    context.occurrenceWords.clear();
    for (size_t length = 16; length >= 2; length /= 2)
    {
        const size_t first = context.rootOccurrences.size();
        for (size_t copy = 0; copy < 2; ++copy)
        {
            EdgeMask mask;
            for (size_t bit = copy * length; bit < (copy + 1) * length; ++bit)
                mask.set(bit);
            context.rootOccurrences.push_back({context.occurrenceWords.size(), 0});
            for (size_t word = 0; word < EdgeMask::activeWordCount(); ++word)
                context.occurrenceWords.push_back(mask.activeWord(word));
        }
        context.rootJobs.push_back({first, first + 1, 0,
            static_cast<uint32_t>(length)});
    }
    const size_t originalJobs = context.rootJobs.size();
    {
        WorkerContext worker(context);
        const auto disabled = runParallelGreedyBootstrap(context, 0, worker, 0);
        assert(disabled.jobsVisited == 0 && disabled.steps == 0);
        const auto noJobs = runParallelGreedyBootstrap(context, 0, worker, 32, 0);
        assert(noJobs.jobsVisited == 0 && noJobs.steps == 0);
        const auto noChecks = runParallelGreedyBootstrap(context, 0, worker, 32, 32, 0);
        assert(noChecks.jobsVisited == 0 && noChecks.steps == 0);
        const auto pastEnd = runParallelGreedyBootstrap(
            context, context.rootJobs.size(), worker
        );
        assert(pastEnd.jobsVisited == 0 && pastEnd.steps == 0);
        assert(worker.assemblyIndex == context.rootAssemblyIndex);
    }
    {
        WorkerContext worker(context);
        const auto one = runParallelGreedyBootstrap(context, 0, worker, 1);
        assert(one.jobsVisited == 1 && one.steps == 1);
        assert(worker.assemblyIndex < context.rootAssemblyIndex);
        assert(worker.root.sumDupBonds == 0);
        assert(worker.root.fragments.size() == 1);
        assert(worker.root.fragments[0].mask.count() == context.bondCount);
    }
    {
        WorkerContext worker(context);
        const auto oneJob = runParallelGreedyBootstrap(context, 0, worker, 32, 1);
        assert(oneJob.jobsVisited == 1 && oneJob.steps == 1);
    }
    {
        WorkerContext worker(context);
        const auto bounded = runParallelGreedyBootstrap(context, 0, worker, 32, 3, 2);
        assert(bounded.jobsVisited <= 3 && bounded.fragmentChecks <= 2);
    }
    {
        WorkerContext worker(context);
        const auto greedy = runParallelGreedyBootstrap(context, 0, worker);
        assert(greedy.steps > 1 && greedy.steps <= 32);
        assert(greedy.jobsVisited <= 65536 && greedy.fragmentChecks <= 262144);
        assert(context.rootJobs.size() == originalJobs);
        assert(worker.root.sumDupBonds == 0);
        assert(worker.root.assemblyIndex() == context.rootAssemblyIndex);
        // The greedy candidate remains a feasible partition of disjoint masks.
        EdgeMask seen;
        for (const assemblyFragment &fragment : worker.candidate.fragments)
        {
            assert((seen & fragment.mask).none());
            assert(fragment.edgeCount == static_cast<int>(fragment.mask.count()));
            seen |= fragment.mask;
        }
        assert(worker.assemblyIndex == worker.candidate.assemblyIndex());
        const int best = worker.assemblyIndex;
        const auto repeated = runParallelGreedyBootstrap(context, 0, worker);
        assert(repeated.steps == greedy.steps && worker.assemblyIndex == best);
        searchCancellationFlag.store(true);
        const auto cancelled = runParallelGreedyBootstrap(context, 0, worker);
        searchCancellationFlag.store(false);
        assert(cancelled.jobsVisited == 0 && cancelled.steps == 0);
        assert(worker.candidate.assemblyIndex() == best);
    }
    clearParallelWorkerMasks();
    return 0;
}
