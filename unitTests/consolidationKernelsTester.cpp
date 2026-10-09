// Exercise shared kernels through their public callers in optimized builds too.
#ifdef NDEBUG
#undef NDEBUG
#endif
#define PARALLELASSEMBLYCPP_NO_MAIN
#include "../src/main.cpp"

#include <cassert>

namespace
{
void configureKernelMaskDomain(size_t width)
{
    bitsetHashTable.clear();
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(width);
    std::construct_at(std::addressof(allEdges));
}

void testVirtualChildBounds()
{
    for (int selected = 2; selected <= 12; ++selected)
    for (int firstEdges = selected; firstEdges <= 30; ++firstEdges)
    for (int secondEdges = selected; secondEdges <= 30; ++secondEdges)
    for (const bool sameParent : {false, true})
    {
        if (sameParent && firstEdges < 2 * selected) continue;
        assemblyState parent;
        parent.appendFragment(EdgeMask{}, firstEdges);
        parent.appendFragment(EdgeMask{}, secondEdges);
        parent.appendFragment(EdgeMask{}, 7);
        IntegerVector totals;
        buildUnrestrictedDupBondTotals(parent, max(2, selected - 1), totals);
        const int second = sameParent ? 0 : 1;
        validMatchings matching(EdgeMaskView{}, EdgeMaskView{}, 0, second, selected);

        // Independently construct the virtual child's edge counts, then apply
        // the ceiling-division bound directly instead of using parent deltas.
        const array<int, 4> childEdges{
            selected,
            firstEdges - (sameParent ? 2 : 1) * selected,
            secondEdges - (sameParent ? 0 : selected),
            7
        };
        int expected = numeric_limits<int>::min();
        for (int size = 2; size <= max(2, selected - 1); ++size)
        {
            int savings = -assembly_bounds::scalarLowerBound(size);
            for (const int edges : childEdges)
                savings += edges - (edges + size - 1) / size;
            expected = max(expected, savings);
        }
        assert(pairSpecificGenericBound(parent, matching, totals) == expected);
        assert(pairSpecificGenericBound(parent, selected, 0, second, totals) == expected);
        assert(pairSpecificGenericBound(parent, selected, second, 0, totals) == expected);
    }
}

void testTouchedAtomTraversal()
{
    for (const size_t width : {63U, 64U, 65U, 127U, 128U, 129U, 513U})
    {
        configureKernelMaskDomain(width);
        ufdsSplit split;
        split.elements.resize(643);
        split.reset();
        // Discover overflow atom words out of order. Components still emerge
        // in atom order in every specialized edge-mask representation.
        split.doubleInsert(641, 642, static_cast<int>(width - 1));
        split.insert(640, 641, static_cast<int>(width - 2));
        split.doubleInsert(512, 513, 5);
        split.insert(511, 512, 6);
        split.doubleInsert(64, 65, 4); // One edge: deliberately omitted.
        split.doubleInsert(0, 1, 0);
        split.insert(2, 1, 1);
        split.merge(0, 2, 2); // Cycle edge comes from extraVals.
        vector<assemblyFragment> fragments;
        split.splitWithBuffers(fragments);
        assert(fragments.size() == 3);
        const array<vector<size_t>, 3> expected{{{0, 1, 2}, {5, 6}, {width - 2, width - 1}}};
        for (size_t component = 0; component < expected.size(); ++component)
        {
            assert(fragments[component].edgeCount == expected[component].size());
            assert(fragments[component].mask.count() == expected[component].size());
            assert(fragments[component].connected);
            assert(fragments[component].canonicalId == unknownCanonicalId);
            for (const size_t edge : expected[component])
                assert(fragments[component].mask.test(edge));
        }
        // Exercise reset and the generation-wrap branch with stale high words.
        fragments.clear();
        split.generation = numeric_limits<uint32_t>::max();
        split.reset();
        split.doubleInsert(512, 513, 7);
        split.insert(511, 512, 8);
        split.splitWithBuffers(fragments);
        assert(fragments.size() == 1);
        assert(fragments.front().mask.count() == 2);
        assert(fragments.front().mask.test(7));
        assert(fragments.front().mask.test(8));
    }
}

void testMatchingTraversalAndPairability()
{
    configureKernelMaskDomain(96);
    vector<potentialDuplicate> occurrences;
    for (const auto &edges : {array<size_t, 2>{0, 64}, {1, 65}, {0, 66}, {0, 67}})
    {
        EdgeMask mask;
        for (const size_t edge : edges) mask.set(edge);
        occurrences.emplace_back(std::move(mask), 0, static_cast<int>(occurrences.size()));
    }
    duplicateSet<potentialDuplicate> duplicates;
    duplicates.bind(2, 1, occurrences, true, false);
    vector<pair<size_t, size_t>> expected;
    for (size_t first = occurrences.size() - 1; first > 0;)
    {
        --first;
        for (size_t second = occurrences.size(); second > first + 1;)
        {
            --second;
            bool overlaps = false;
            for (size_t edge = 0; edge < EdgeMask::size(); ++edge)
                overlaps |= occurrences[first].mask[edge] && occurrences[second].mask[edge];
            if (!overlaps) expected.emplace_back(first, second);
        }
    }
    vector<pair<size_t, size_t>> filtered;
    assert(duplicates.visitMatchingsInReverse(
        [&](validMatchings &, size_t first, size_t second)
        {
            filtered.emplace_back(first, second);
            return false;
        },
        [](validMatchings &) { return true; }
    ));
    assert(filtered == expected);
    vector<pair<EdgeMask, EdgeMask>> unfiltered;
    assert(duplicates.visitMatchingsInReverse([&](validMatchings &matching)
    {
        unfiltered.emplace_back(matching.first.toMask(), matching.second.toMask());
        return true;
    }));
    assert(unfiltered.size() == expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        assert(unfiltered[i].first == occurrences[expected[i].first].mask);
        assert(unfiltered[i].second == occurrences[expected[i].second].mask);
    }
    size_t visits = 0;
    assert(!duplicates.visitMatchingsInReverse([&](validMatchings &)
    {
        return ++visits < 2;
    }));
    assert(visits == 2);
    visits = 0;
    assert(duplicates.visitMatchingsInReverse(
        [](validMatchings &, size_t, size_t) { return true; },
        [&](validMatchings &) { ++visits; return true; }
    ));
    assert(visits == 0);

    // Append an occurrence overlapping every prior one so the callback also
    // has to observe a dead occurrence, after the living ones in input order.
    EdgeMask overlap;
    for (size_t edge = 0; edge < 2; ++edge) overlap.set(edge);
    occurrences.emplace_back(std::move(overlap), 0, 4);
    duplicates.bind(2, 1, occurrences, true, false);
    vector<uint8_t> scratch;
    vector<bool> alive;
    assert(duplicates.visitOccurrencePairability(scratch,
        [&](potentialDuplicate &occurrence, bool pairable)
        {
            assert(occurrence.duplicateIndex == static_cast<int>(alive.size()));
            alive.push_back(pairable);
            return true;
        }
    ));
    assert(alive == vector<bool>({true, true, true, true, false}));
    visits = 0;
    assert(!duplicates.visitOccurrencePairability(scratch,
        [&](potentialDuplicate &, bool) { return ++visits < 2; }
    ));
    assert(visits == 2);

    occurrences.back().fragmentIndex = 1;
    duplicates.bind(2, 2, occurrences, true, true);
    scratch.assign(1, 123);
    alive.clear();
    assert(duplicates.visitOccurrencePairability(scratch,
        [&](potentialDuplicate &, bool pairable) { alive.push_back(pairable); return true; }
    ));
    assert(alive == vector<bool>(occurrences.size(), true));
    assert(scratch == vector<uint8_t>{123}); // Multi-fragment shortcut.

    searchCancellationFlag.store(true);
    visits = 0;
    assert(!duplicates.visitMatchingsInReverse([&](validMatchings &) { ++visits; return true; }));
    assert(!duplicates.visitOccurrencePairability(scratch,
        [&](potentialDuplicate &, bool) { ++visits; return true; }
    ));
    assert(visits == 0);
    searchCancellationFlag.store(false);
}

void testMolecularInputInitialization()
{
    for (const size_t width : {0U, 2U, 65U, 129U, 64U, 0U})
    {
        molGraph graph;
        for (size_t vertex = 0; vertex <= width; ++vertex) graph.addAtom("C");
        for (size_t edge = 0; edge < width; ++edge)
            graph.addBond(static_cast<int>(edge), static_cast<int>(edge + 1), 1);
        const bool hasUniqueBond = width >= 2;
        if (hasUniqueBond)
        {
            graph.addAtom("O");
            graph.addBond(0, static_cast<int>(width + 1), 2);
        }
        graph.addAtom("N"); // Isolated atoms remain in the input convention.
        vector<MoleculeEdge> removed;
        runtimeLimitReached = true;
        enumerationLimitReached = true;
        intermediateAssemblyIndices.emplace_back(1, 123);
        prepareMolecularSearchInput(graph, removed, clock_t{17});
        assert(startTime == clock_t{17});
        assert(!runtimeLimitReached && !enumerationLimitReached);
        assert(intermediateAssemblyIndices.empty());
        assert(bitsetHashTable.empty() && graphHashMap.empty());
        assert(sharedTargetMolecule == nullptr && sharedUniverseEdgeList == nullptr);
        assert(sharedCanonicalRegistry == nullptr && sharedAssemblyStates == nullptr);
        assert(originalMolecule.atoms.size() == graph.atoms.size());
        assert(originalEdgeList.size() == width + hasUniqueBond);
        assert(totalBonds == width + hasUniqueBond);
        assert(disjointFragments == 2);
        assert(removed.size() == static_cast<size_t>(hasUniqueBond));
        assert(universeEdgeList.size() == width);
        assert(EdgeMask::size() == width);
        assert(AtomMask::size() == graph.atoms.size());
        assert(allEdges.none());
        allEdges.set();
        assert(allEdges.count() == width);
        if (width != 0)
        {
            assert(canonise(allEdges) >= 0);
            assert(!bitsetHashTable.empty());
        }
    }
}

void testBondedComponentCounts()
{
    for (int isolated = 0; isolated < 4; ++isolated)
    for (int bonded = 0; bonded < 4; ++bonded)
    {
        molGraph graph;
        for (int atomIndex = 0; atomIndex < isolated + 3 * bonded; ++atomIndex)
            graph.addAtom("C");
        for (int component = 0; component < bonded; ++component)
        {
            const int start = isolated + 3 * component;
            graph.addBond(start, start + 1, 1);
            graph.addBond(start + 1, start + 2, 2);
        }
        assert(graph.disjointFragments() == isolated + bonded);
        assert(graphRepair::implementation::componentCount(graph) == bonded);
    }
}
} // namespace

int main()
{
    suppressSearchOutput = true;
    maximumRuntimeTicks = numeric_limits<unsigned long long>::max();
    testVirtualChildBounds();
    testTouchedAtomTraversal();
    testMatchingTraversalAndPairability();
    testBondedComponentCounts();
    testMolecularInputInitialization();
    configureKernelMaskDomain(64);
}
