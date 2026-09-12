#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <csignal>
#include <ctime>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#define PARALLELASSEMBLYCPP_SEARCH_LOCAL thread_local
using namespace std;
using IntegerVector = vector<int>;
using BooleanVector = vector<bool>;
using IntegerPair = pair<int, int>;

#include "../src/activeWordMask.h"
#include "../src/globalPrimitives.h"
#include "../src/searchTelemetry.h"
#include "../src/matchingBoundRefresh.h"

namespace
{

void resetCounters()
{
    searchTelemetryEnabled = true;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    searchTelemetry = SearchTelemetryState{};
#endif
}

void expectPolls(uint64_t polls, uint64_t refreshes)
{
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    assert(searchTelemetry.counters.matchingBoundRefreshPolls == polls);
    assert(searchTelemetry.counters.matchingBoundRefreshes == refreshes);
#else
    static_cast<void>(polls);
    static_cast<void>(refreshes);
#endif
}

template<bool enabled>
unsigned int runLongScan()
{
    resetCounters();
    std::atomic<int> sharedBest{20};
    int best = 20;
    matchingBoundRefresh<enabled> refresh(&sharedBest);
    assert(refresh.atClassBoundary(best) == 20);

    unsigned int candidatesExpanded = 0;
    // A worker finds a better incumbent while this worker stays in the same
    // class, without reaching another recursive entry or class boundary.
    for (unsigned int i = 0; i < 4 * refresh.interval; ++i)
    {
        if (i == refresh.interval / 2)
            sharedBest.store(10, std::memory_order_relaxed);
        const int previousBest = refresh.beforeMatchingBound(best);
        constexpr int lowerBound = 12;
        if (lowerBound >= best)
        {
            refresh.recordPrune(
                lowerBound, previousBest, best, matchingBoundWork::candidate
            );
            continue;
        }
        ++candidatesExpanded;
    }

    if constexpr (enabled)
    {
        assert(best == 10);
        expectPolls(5, 1);
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        // Only the prune immediately following the refresh is attributable;
        // later avoided expansions must not inflate this measurement.
        assert(searchTelemetry.counters.matchingBoundCandidatesPruned == 1);
#endif
    }
    else
    {
        assert(best == 20);
        expectPolls(0, 0);
    }
    return candidatesExpanded;
}

void testLongScan()
{
    const unsigned int withoutRefresh = runLongScan<false>();
    const unsigned int withRefresh = runLongScan<true>();
    assert(withoutRefresh == 4 * matchingBoundRefresh<true>::interval);
    assert(withRefresh == matchingBoundRefresh<true>::interval - 1);
    assert(withRefresh < withoutRefresh);
}

void testClassRefreshAndCadenceReset()
{
    resetCounters();
    std::atomic<int> sharedBest{18};
    int best = 20;
    matchingBoundRefresh<true> refresh(&sharedBest);
    assert(refresh.atClassBoundary(best) == 20);
    assert(best == 18);
    for (unsigned int i = 0; i < refresh.interval / 2; ++i)
        assert(refresh.beforeMatchingBound(best) == 18);
    expectPolls(1, 1);

    // A new class refreshes immediately and starts a fresh full interval.
    sharedBest.store(16, std::memory_order_relaxed);
    assert(refresh.atClassBoundary(best) == 18);
    assert(best == 16);
    sharedBest.store(14, std::memory_order_relaxed);
    for (unsigned int i = 1; i < refresh.interval; ++i)
    {
        assert(refresh.beforeMatchingBound(best) == 16);
        assert(best == 16);
    }
    expectPolls(2, 2);
    assert(refresh.beforeMatchingBound(best) == 16);
    assert(best == 14);
    expectPolls(3, 3);
}

void testRefreshNeverWeakensBound()
{
    resetCounters();
    std::atomic<int> sharedBest{20};
    int best = 20;
    matchingBoundRefresh<true> refresh(&sharedBest);
    assert(refresh.atClassBoundary(best) == 20);
    assert(best == 20);
    sharedBest.store(22, std::memory_order_relaxed);
    assert(refresh.atClassBoundary(best) == 20);
    assert(best == 20);
    expectPolls(2, 0);

    // A local improvement can also precede publication to the shared bound.
    best = 10;
    for (unsigned int i = 0; i < refresh.interval; ++i)
        assert(refresh.beforeMatchingBound(best) == 10);
    assert(best == 10);
    expectPolls(3, 0);
    sharedBest.store(8, std::memory_order_relaxed);
    assert(refresh.atClassBoundary(best) == 10);
    assert(best == 8);
    expectPolls(4, 1);
}

template<bool enabled>
void checkInactiveRefresh(std::atomic<int> *sharedBest)
{
    resetCounters();
    int best = 20;
    matchingBoundRefresh<enabled> refresh(sharedBest);
    for (unsigned int classIndex = 0; classIndex < 3; ++classIndex)
    {
        assert(refresh.atClassBoundary(best) == 20);
        for (unsigned int i = 0; i < 2 * refresh.interval; ++i)
        {
            assert(refresh.beforeMatchingBound(best) == 20);
            assert(best == 20);
        }
    }
    expectPolls(0, 0);
}

void testInactiveRefresh()
{
    std::atomic<int> sharedBest{10};
    checkInactiveRefresh<false>(&sharedBest);
    checkInactiveRefresh<false>(nullptr);
    checkInactiveRefresh<true>(nullptr);
}

void testPruneAttribution()
{
    resetCounters();
    constexpr std::array workKinds{
        matchingBoundWork::duplicateClass,
        matchingBoundWork::occurrencePair,
        matchingBoundWork::fragmentPairBlock,
        matchingBoundWork::candidate
    };
    for (const matchingBoundWork work : workKinds)
    {
        // Equal-to-new and intermediate bounds now prune; below-new bounds
        // still need work, and equal-to-old bounds were already prunable.
        matchingBoundRefresh<true>::recordPrune(10, 20, 10, work);
        matchingBoundRefresh<true>::recordPrune(15, 20, 10, work);
        matchingBoundRefresh<true>::recordPrune(9, 20, 10, work);
        matchingBoundRefresh<true>::recordPrune(20, 20, 10, work);
        matchingBoundRefresh<true>::recordPrune(21, 20, 10, work);
        matchingBoundRefresh<true>::recordPrune(15, 10, 10, work);
        matchingBoundRefresh<true>::recordPrune(15, 10, 20, work);
        matchingBoundRefresh<false>::recordPrune(15, 20, 10, work);
    }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    const SearchTelemetryCounters &counters = searchTelemetry.counters;
    assert(counters.matchingBoundClassesPruned == 2);
    assert(counters.matchingBoundPairsPruned == 2);
    assert(counters.matchingBoundBlocksPruned == 2);
    assert(counters.matchingBoundCandidatesPruned == 2);
#endif
    expectPolls(0, 0);
}

void testRuntimeTelemetryDisabled()
{
    resetCounters();
    searchTelemetryEnabled = false;
    std::atomic<int> sharedBest{10};
    int best = 20;
    matchingBoundRefresh<true> refresh(&sharedBest);
    assert(refresh.atClassBoundary(best) == 20);
    assert(best == 10);
    sharedBest.store(8, std::memory_order_relaxed);
    for (unsigned int i = 0; i < refresh.interval; ++i)
        refresh.beforeMatchingBound(best);
    assert(best == 8);
    refresh.recordPrune(8, 10, 8, matchingBoundWork::candidate);
    expectPolls(0, 0);
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    assert(searchTelemetry.counters.matchingBoundCandidatesPruned == 0);
#endif
    searchTelemetryEnabled = true;
}

#ifdef ASSEMBLY_ENABLE_TELEMETRY
void testTelemetryAggregationAndJson()
{
    SearchTelemetryCounters aggregate;
    SearchTelemetryCounters worker;
    worker.matchingBoundRefreshPolls = 2;
    worker.matchingBoundRefreshes = 3;
    worker.matchingBoundClassesPruned = 5;
    worker.matchingBoundPairsPruned = 7;
    worker.matchingBoundBlocksPruned = 11;
    worker.matchingBoundCandidatesPruned = 13;
    addSearchTelemetryCounters(aggregate, worker);
    addSearchTelemetryCounters(aggregate, worker);
    assert(aggregate.matchingBoundRefreshPolls == 4);
    assert(aggregate.matchingBoundRefreshes == 6);
    assert(aggregate.matchingBoundClassesPruned == 10);
    assert(aggregate.matchingBoundPairsPruned == 14);
    assert(aggregate.matchingBoundBlocksPruned == 22);
    assert(aggregate.matchingBoundCandidatesPruned == 26);

    std::ostringstream output;
    writeAllSearchTelemetryCounters(output, aggregate, "");
    const std::string json = output.str();
    for (const char *field : {
        "\"matching_bound_refresh_polls\": 4,",
        "\"matching_bound_refreshes\": 6,",
        "\"matching_bound_classes_pruned\": 10,",
        "\"matching_bound_pairs_pruned\": 14,",
        "\"matching_bound_blocks_pruned\": 22,",
        "\"matching_bound_candidates_pruned\": 26,"
    })
    {
        assert(json.find(field) != std::string::npos);
    }
}
#endif

}

int main()
{
    testLongScan();
    testClassRefreshAndCadenceReset();
    testRefreshNeverWeakensBound();
    testInactiveRefresh();
    testPruneAttribution();
    testRuntimeTelemetryDisabled();
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    testTelemetryAggregationAndJson();
#endif
}
