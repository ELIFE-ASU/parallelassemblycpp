// Keep oracle checks active in Release/CI builds as well.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../src/additionChainBounds.h"

#include <cassert>
#include <climits>
#include <iostream>
#include <set>

namespace
{

// This scalar oracle searches complete increasing chains, independently of the
// production lookup table, and allows every pair of earlier elements.
bool scalarReachable(int target, int remaining, std::vector<int> &chain)
{
    if (chain.back() == target) return true;
    if (remaining == 0
        || (static_cast<std::uint64_t>(chain.back()) << remaining)
            < static_cast<std::uint64_t>(target))
        return false;
    std::set<int, std::greater<int>> nextValues;
    for (std::size_t i = 0; i < chain.size(); ++i)
        for (std::size_t j = 0; j <= i; ++j)
        {
            const int value = chain[i] + chain[j];
            if (value > chain.back() && value <= target) nextValues.insert(value);
        }
    for (const int value : nextValues)
    {
        chain.push_back(value);
        const bool found = scalarReachable(target, remaining - 1, chain);
        chain.pop_back();
        if (found) return true;
    }
    return false;
}

void testScalarTable()
{
    assert(assembly_bounds::scalarLowerBound(0) == 0);
    for (int target = 1; target <= 256; ++target)
    {
        std::vector<int> chain{1};
        int exact = 0;
        while ((1 << exact) < target) ++exact;
        while (!scalarReachable(target, exact, chain)) ++exact;
        assert(assembly_bounds::scalarLowerBound(target) == exact);
    }
    assert(assembly_bounds::scalarLowerBound(7) == 4);
    assert(assembly_bounds::scalarLowerBound(15) == 5);
    assert(assembly_bounds::scalarLowerBound(255) == 10);
    assert(assembly_bounds::scalarLowerBound(256) == 8);
    assert(assembly_bounds::scalarLowerBound(257) == 9);
    assert(assembly_bounds::scalarLowerBound(1024) == 10);
    assert(assembly_bounds::scalarLowerBound(INT_MAX) == 31);
}

using Vector = std::vector<int>;
using State = std::vector<Vector>;

// Independent breadth-first oracle over sets of available vectors. Unlike the
// production IDDFS, this neither encodes coordinates nor orders generated
// vectors. Merging equivalent sets still enumerates every possible pair sum.
int vectorOracle(const Vector &target)
{
    State initial;
    for (std::size_t c = 0; c < target.size(); ++c)
    {
        Vector unit(target.size());
        unit[c] = 1;
        initial.push_back(unit);
        if (unit == target) return 0;
    }
    std::sort(initial.begin(), initial.end());
    std::set<State> frontier{initial};
    for (int depth = 1;; ++depth)
    {
        std::set<State> next;
        for (const State &state : frontier)
            for (std::size_t i = 0; i < state.size(); ++i)
                for (std::size_t j = 0; j <= i; ++j)
                {
                    Vector sum(target.size());
                    bool fits = true;
                    for (std::size_t c = 0; c < target.size(); ++c)
                    {
                        sum[c] = state[i][c] + state[j][c];
                        if (sum[c] > target[c]) fits = false;
                    }
                    if (!fits) continue;
                    if (sum == target) return depth;
                    if (std::binary_search(state.begin(), state.end(), sum)) continue;
                    State successor = state;
                    successor.insert(std::lower_bound(
                        successor.begin(), successor.end(), sum
                    ), sum);
                    next.insert(std::move(successor));
                }
        assert(!next.empty());
        frontier = std::move(next);
    }
}

void testVectorOracle()
{
    // All partitions of every total through 8 using at most three coordinates.
    // Every ordering of the named coordinates has the same optimum.
    assembly_bounds::VectorBoundCache generous(256, 1000000);
    for (int total = 1; total <= 8; ++total)
        for (int a = 1; a <= total; ++a)
            for (int b = 0; b <= a; ++b)
                for (int c = 0; c <= b; ++c)
                {
                    if (a + b + c != total) continue;
                    Vector counts{a};
                    if (b != 0) counts.push_back(b);
                    if (c != 0) counts.push_back(c);
                    const int exact = vectorOracle(counts);
                    assert(generous.lowerBound(counts) == exact);
                    for (int threshold = 0; threshold <= exact + 1; ++threshold)
                    {
                        assembly_bounds::VectorBoundCache thresholded(0, 1000000);
                        const int certificate = thresholded.lowerBound(counts, threshold);
                        assert(certificate <= exact);
                        assert(certificate >= std::min(exact, threshold));
                    }
                    int previous = -1;
                    for (const std::size_t budget : {0U, 1U, 7U, 31U, 127U, 4000U})
                    {
                        assembly_bounds::VectorBoundCache bounded(0, budget);
                        const int bound = bounded.lowerBound(counts);
                        assert(bound >= assembly_bounds::scalarLowerBound(total));
                        assert(bound >= static_cast<int>(counts.size()) - 1);
                        assert(bound <= exact);
                        assert(bound >= previous);
                        previous = bound;
                    }
                }
}

void testFourCoordinateOracle()
{
    assembly_bounds::VectorBoundCache generous(256, 1000000);
    for (const Vector &counts : std::vector<Vector>{
        {2, 1, 1, 1}, {2, 2, 1, 1}, {3, 1, 1, 1}, {2, 2, 2, 1}, {2, 2, 2, 2}
    })
    {
        const int exact = vectorOracle(counts);
        assert(generous.lowerBound(counts) == exact);
        for (const std::size_t budget : {0U, 1U, 31U, 4000U})
        {
            assembly_bounds::VectorBoundCache bounded(0, budget);
            assert(bounded.lowerBound(counts) <= exact);
        }
    }
}

void testThresholdCacheAndSupportBound()
{
    // Combining five distinct units takes four additions; doubling once gives
    // (2,2,2,2,2). A five-source tree with only four additions cannot repeat a
    // source, independently proving that five is the exact optimum.
    assembly_bounds::VectorBoundCache analyticalOnly(0, 0);
    assert(analyticalOnly.lowerBound({2, 2, 2, 2, 2}) == 5);
    assert(analyticalOnly.lowerBound({1, 1, 1, 1, 2}) == 5);
    assert(analyticalOnly.lowerBound({1, 1, 1, 1, 1}) == 4);
    assembly_bounds::VectorBoundCache thresholded;
    assert(thresholded.lowerBound({2, 2, 2, 2, 2}, 5) == 5);
    assert(thresholded.lowerBound({1, 3}, 2) == 2);
    // Cache entries are certificates, not promises to exhaust a new threshold.
    assert(thresholded.lowerBound({3, 1}) == 2);
    assert(thresholded.lowerBound({3, 1}, 3) == 2);
    thresholded.clear();
    assert(thresholded.lowerBound({3, 1}, 3) == 3);
    assert(thresholded.lowerBound({3, 1}, 2) == 3);
}

void testVectorExamplesAndCanonicalCache()
{
    assembly_bounds::VectorBoundCache cache;
    assert(cache.lowerBound({}) == 0);
    assert(cache.lowerBound({0, 0}) == 0);
    assert(cache.lowerBound({1, 0}) == 0);
    cache.clear();
    assert(cache.size() == 0);
    assert(cache.lowerBound({3, 1}) == 3); // AAAB beats the scalar floor 2.
    assert(cache.size() == 1);
    assert(cache.lowerBound({0, 1, 0, 3, 0}) == 3);
    assert(cache.size() == 1);
    assert(cache.lowerBound({2, 2}) == 2);
    assert(cache.lowerBound({3, 3}) == 3);
    assert(cache.lowerBound({1, 1, 1, 1, 1, 1}) == 5);
    assert(cache.lowerBound({16}) == 4);

    // Large vectors use certified projections without exact search. The total
    // can exceed int even though each count is a valid int.
    assert(cache.lowerBound({INT_MAX, INT_MAX}) == 32);
    assert(cache.lowerBound(Vector(300, 1)) == 299);
    assert(cache.lowerBound({255, 1}) == 10); // total=256 would only give 8.

    bool rejected = false;
    try { static_cast<void>(cache.lowerBound({2, -1})); }
    catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);

    assembly_bounds::VectorBoundCache tiny(2);
    for (int i = 1; i < 20; ++i)
    {
        assert(tiny.lowerBound({i}) == assembly_bounds::scalarLowerBound(i));
        assert(tiny.size() <= 2);
    }
    assert(tiny.lowerBound({3, 1}) == 3);
    assembly_bounds::VectorBoundCache noCache(0);
    assert(noCache.lowerBound({3, 1}) == 3);
    assert(noCache.size() == 0);
    // Independent instances must not share mutable cache state.
    assembly_bounds::VectorBoundCache second;
    assert(second.size() == 0);
    assert(second.lowerBound({3, 1}) == 3);
    cache.clear();
    assert(second.size() == 1);
}

} // namespace

int main()
{
    testScalarTable();
    testVectorOracle();
    testFourCoordinateOracle();
    testThresholdCacheAndSupportBound();
    testVectorExamplesAndCanonicalCache();
    std::cout << "Addition-chain bound tests passed\n";
}
