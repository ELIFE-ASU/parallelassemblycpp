#pragma once

enum class matchingBoundWork
{
    duplicateClass,
    occurrencePair,
    fragmentPairBlock,
    candidate
};

/** Poll only the incumbent value: it does not publish or acquire a pathway. */
template<bool enabled>
class matchingBoundRefresh
{
    std::atomic<int> *shared;
    unsigned int remaining = interval;

    int refresh(int &bestAssemblyIndex)
    {
        const int previous = bestAssemblyIndex;
        if constexpr (enabled)
        {
            if (shared != nullptr)
            {
                const int observed = shared->load(std::memory_order_relaxed);
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                if (searchTelemetryEnabled) [[unlikely]]
                    ++searchTelemetry.counters.matchingBoundRefreshPolls;
#endif
                if (observed < bestAssemblyIndex)
                {
                    bestAssemblyIndex = observed;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                    if (searchTelemetryEnabled) [[unlikely]]
                        ++searchTelemetry.counters.matchingBoundRefreshes;
#endif
                }
            }
        }
        return previous;
    }

public:
    // Bound evaluations, rather than recursive entries, advance this cadence.
    static constexpr unsigned int interval = 64;

    explicit matchingBoundRefresh(std::atomic<int> *incumbent) : shared(incumbent) {}

    int atClassBoundary(int &bestAssemblyIndex)
    {
        remaining = interval;
        return refresh(bestAssemblyIndex);
    }

    int beforeMatchingBound(int &bestAssemblyIndex)
    {
        if constexpr (enabled)
        {
            if (shared != nullptr && --remaining == 0)
            {
                remaining = interval;
                return refresh(bestAssemblyIndex);
            }
        }
        return bestAssemblyIndex;
    }

    // Count only immediately demonstrated savings. Later prunes and changes
    // in the search tree cannot be attributed to a refresh without replaying it.
    static void recordPrune(
        int lowerBound,
        int previousBest,
        int currentBest,
        matchingBoundWork work
    )
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if constexpr (enabled)
        {
            if (
                searchTelemetryEnabled &&
                lowerBound >= currentBest && lowerBound < previousBest
            ) [[unlikely]]
            {
                switch (work)
                {
                    case matchingBoundWork::duplicateClass:
                        ++searchTelemetry.counters.matchingBoundClassesPruned;
                        break;
                    case matchingBoundWork::occurrencePair:
                        ++searchTelemetry.counters.matchingBoundPairsPruned;
                        break;
                    case matchingBoundWork::fragmentPairBlock:
                        ++searchTelemetry.counters.matchingBoundBlocksPruned;
                        break;
                    case matchingBoundWork::candidate:
                        ++searchTelemetry.counters.matchingBoundCandidatesPruned;
                        break;
                }
            }
        }
#else
        static_cast<void>(lowerBound);
        static_cast<void>(previousBest);
        static_cast<void>(currentBest);
        static_cast<void>(work);
#endif
    }
};
