#pragma once

struct parallelGreedyBootstrapResult
{
    size_t jobsVisited = 0;
    size_t fragmentChecks = 0;
    size_t steps = 0;
};

[[nodiscard]] inline bool parallelGreedyBootstrapEnabled()
{
    const char *setting = std::getenv("PARALLELASSEMBLYCPP_GREEDY_BOOTSTRAP");
    if (setting == nullptr || std::string_view(setting) == "0") return false;
    if (std::string_view(setting) == "1") return true;
    throw std::invalid_argument(
        "PARALLELASSEMBLYCPP_GREEDY_BOOTSTRAP must be 0 or 1"
    );
}

/**
 * Follow one feasible largest-first path through the immutable root pairs.
 * Fragmentation only splits/removes masks, so a rejected pair cannot become
 * usable later. Re-find parents in the current fragments before each step.
 * Hard job, step, and parent-scan budgets bound this experiment independently
 * of the root frontier size. It neither consumes jobs nor caches searched
 * states; exhaustive search still starts from the original root.
 */
parallelGreedyBootstrapResult runParallelGreedyBootstrap(
    const SearchContext &context,
    size_t jobIndex,
    WorkerContext &worker,
    size_t maximumSteps = 32,
    size_t maximumJobs = 65536,
    size_t maximumFragmentChecks = 262144
)
{
    parallelGreedyBootstrapResult result;
    assemblyState next;
    assemblyState *current = &worker.root;
    while (
        jobIndex < context.rootJobs.size() &&
        result.jobsVisited < maximumJobs && result.steps < maximumSteps &&
        result.fragmentChecks < maximumFragmentChecks
    )
    {
        if (searchShouldStop()) break;
        const rootJobDescriptor &job = context.rootJobs[jobIndex++];
        ++result.jobsVisited;
        EdgeMask first = reconstructRootOccurrence(context, job.firstOccurrence);
        EdgeMask second = reconstructRootOccurrence(context, job.secondOccurrence);
        int firstParent = -1;
        int secondParent = -1;
        for (size_t index = 0; index < current->fragments.size(); ++index)
        {
            if (result.fragmentChecks == maximumFragmentChecks) return result;
            ++result.fragmentChecks;
            const assemblyFragment &fragment = current->fragments[index];
            if (fragment.edgeCount < static_cast<int>(job.duplicateSize)) continue;
            if (fragment.mask.contains(first)) firstParent = static_cast<int>(index);
            if (fragment.mask.contains(second)) secondParent = static_cast<int>(index);
            if (firstParent >= 0 && secondParent >= 0) break;
        }
        if (firstParent < 0 || secondParent < 0) continue;
        validMatchings matching(
            first, second, firstParent, secondParent,
            static_cast<int>(job.duplicateSize)
        );
        next.clearFragments();
        fragmentAssemblyStateWithoutCanonisationWithWorkspace(
            *current, matching, job.canonicalId, next, worker.fragmentation
        );
        if (searchShouldStop()) break;
        next.sumDupBonds = current->sumDupBonds + matching.maximumFragmentSize - 1;
        swap(worker.candidate, next);
        current = &worker.candidate;
        ++result.steps;
        recordImprovedAssemblyIndex<false>(
            worker.candidate, worker.assemblyIndex, worker.search
        );
    }
    return result;
}
