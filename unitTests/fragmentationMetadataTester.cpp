#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <csignal>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;
using vi = vector<int>;
using vb = vector<bool>;
using pii = pair<int, int>;

#include "../v5/activeWordMask.h"

constexpr int ceilLog2(int value)
{
    return std::bit_width(static_cast<unsigned int>(value - 1));
}

#include "../v5/globalPrimitives.h"
#include "../v5/ufds.h"
#include "../v5/molGraph.h"
#include "../v5/treeCanon.h"
#include "../v5/cyclicCanon.h"
#include "../v5/assemblyState.h"
#include "../v5/graphHashes.h"
#include "../v5/dagEnumeration.h"
#include "../v5/duplicateMatching.h"
#include "../v5/fragmentation.h"

void require(bool condition)
{
    if (!condition) abort();
}

#undef assert
#define assert(condition) require(static_cast<bool>(condition))

EdgeMask makeMask(initializer_list<size_t> edges)
{
    EdgeMask result;
    for (const size_t edge : edges) result.set(edge);
    return result;
}

void configurePathGraph(size_t edgeCount)
{
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    univEdgeList.clear();
    targetMolecule = molGraph();
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(edgeCount);
    std::construct_at(std::addressof(allEdges));
    string atomType = "C";
    for (size_t atomIndex = 0; atomIndex <= edgeCount; atomIndex++)
        targetMolecule.addAtom(atomType);
    for (size_t edge = 0; edge < edgeCount; edge++)
        targetMolecule.addBond(edge, edge + 1, 1);
    univEdgeList = targetMolecule.writeEdgeList();
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
    AtomMask::configure(edgeCount + 1);
    totalBonds = edgeCount;
}

vector<int> resolveFirstCacheReplay(
    const EdgeMask &mask,
    ufdsMaskWorkspace &workspace
)
{
    assemblyState output;
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(mask, output.fragments, workspace);
    assert(!output.fragments.empty());

    output.clearFragments();
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(mask, output.fragments, workspace);
    assert(!output.fragments.empty());

    output.clearFragments();
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(mask, output.fragments, workspace);
    assert(!output.fragments.empty());
    for (const assemblyFragment &fragment : output.fragments)
        assert(fragment.canonicalId == unknownCanonicalId);

    vi key;
    assert(canoniseAssemblyStateAndBuildKey(output, key, workspace));
    vector<int> canonicalIds;
    for (const assemblyFragment &fragment : output.fragments)
        canonicalIds.push_back(fragment.canonicalId);
    return canonicalIds;
}

void admitUnresolvedResidual(
    const EdgeMask &mask,
    ufdsMaskWorkspace &workspace
)
{
    assemblyState output;
    for (int request = 0; request < 2; request++)
    {
        workspace.beginFragmentation();
        output.clearFragments();
        ufdsMaskConstructWithWorkspace(mask, output.fragments, workspace);
        assert(!output.fragments.empty());
    }
}

void testMaskTrimmingMetadata()
{
    assemblyFragment unchanged(makeMask({0, 1, 2}), 3, 123, true);
    assert(!unchanged.retainEdges(makeMask({0, 1, 2, 3})));
    assert(unchanged.mask == makeMask({0, 1, 2}));
    assert(unchanged.edgeCount == 3);
    assert(unchanged.canonicalId == 123);
    assert(unchanged.connected);

    assemblyFragment changed(makeMask({4, 5, 6, 7}), 4, 456, true);
    assert(changed.retainEdges(makeMask({4, 5, 6})));
    assert(changed.mask == makeMask({4, 5, 6}));
    assert(changed.edgeCount == 3);
    assert(changed.canonicalId == unknownCanonicalId);
    assert(!changed.connected);
}

void testFragmentMetadataAndSelectiveCanonisation()
{
    const assemblyFragment maximumId(
        makeMask({30, 31}),
        2,
        numeric_limits<int32_t>::max(),
        true
    );
    assert(maximumId.canonicalId == numeric_limits<int32_t>::max());

    assemblyState target;
    target.appendFragment(makeMask({0, 1, 2, 3, 4, 5}), 6, 50, true);
    target.appendFragment(makeMask({6, 7}), 2, 42, true);
    target.appendFragment(
        makeMask({8, 9, 12, 13}),
        4,
        unknownCanonicalId,
        false
    );

    EdgeMask first = makeMask({0, 1});
    EdgeMask second = makeMask({4, 5});
    validMatchings matching(first, second, 0, 0, 2);
    ufdsMaskWorkspace workspace(33, 32);
    assemblyState child;
    fragmentAssemblyStateWithoutCanonisationWithWorkspace(
        target,
        matching,
        91,
        child,
        workspace
    );

    assert(child.fragments.size() == 5);
    assert(child.fragments[0].mask == first);
    assert(child.fragments[0].canonicalId == 91);
    assert(child.fragments[0].connected);
    assert(child.fragments[1].mask == makeMask({2, 3}));
    assert(child.fragments[1].canonicalId == unknownCanonicalId);
    assert(child.fragments[1].connected);
    assert(child.fragments[2].mask == target.fragments[1].mask);
    assert(child.fragments[2].canonicalId == 42);
    assert(child.fragments[2].connected);
    assert(child.fragments[3].mask == makeMask({8, 9}));
    assert(child.fragments[3].canonicalId == unknownCanonicalId);
    assert(child.fragments[3].connected);
    assert(child.fragments[4].mask == makeMask({12, 13}));
    assert(child.fragments[4].canonicalId == unknownCanonicalId);
    assert(child.fragments[4].connected);

    const size_t initialCachedMaskCount = bitsetHashTable.size();
    vi key;
    assert(canoniseAssemblyStateAndBuildKey(child, key, workspace));
    assert(bitsetHashTable.size() == initialCachedMaskCount + 3);
    assert(key.size() == child.fragments.size());
    assert(key.front() == 91);
    assert(child.fragments[2].canonicalId == 42);
    assert(child.fragments[1].canonicalId >= 0);
    assert(
        child.fragments[1].canonicalId == child.fragments[3].canonicalId
    );
    assert(is_sorted(key.begin() + 1, key.end()));

    const size_t cachedMaskCount = bitsetHashTable.size();
    const size_t canonicalClassCount = graphHashMap.size();
    const vi firstKey = key;
    assert(canoniseAssemblyStateAndBuildKey(child, key, workspace));
    assert(key == firstKey);
    assert(bitsetHashTable.size() == cachedMaskCount);
    assert(graphHashMap.size() == canonicalClassCount);
}

void testMultipleResidualCacheBindings()
{
    ufdsMaskWorkspace workspace(33, 32);
    const EdgeMask identity = makeMask({16, 17, 18, 19});
    const EdgeMask split = makeMask({20, 21, 24, 25});
    const uint64_t identityFingerprint =
        ufdsMaskWorkspace::mixDecompositionFingerprint(
            bitsetLowWordBelow(identity, workspace.edgeCount)
        );
    const uint64_t splitFingerprint =
        ufdsMaskWorkspace::mixDecompositionFingerprint(
            bitsetLowWordBelow(split, workspace.edgeCount)
        );
    assert(
        (identityFingerprint &
            (ufdsMaskWorkspace::decompositionCacheEntryLimit - 1)) !=
        (splitFingerprint &
            (ufdsMaskWorkspace::decompositionCacheEntryLimit - 1))
    );

    admitUnresolvedResidual(identity, workspace);
    admitUnresolvedResidual(split, workspace);

    assemblyState child;
    child.appendFragment(makeMask({0, 1}), 2, 700, true);
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(identity, child.fragments, workspace);
    ufdsMaskConstructWithWorkspace(split, child.fragments, workspace);
    assert(child.fragments.size() == 4);
    for (size_t i = 1; i < child.fragments.size(); i++)
        assert(child.fragments[i].canonicalId == unknownCanonicalId);

    vi key;
    assert(canoniseAssemblyStateAndBuildKey(child, key, workspace));
    vector<int> resolvedIds;
    for (size_t i = 1; i < child.fragments.size(); i++)
        resolvedIds.push_back(child.fragments[i].canonicalId);

    assemblyState replay;
    replay.appendFragment(makeMask({0, 1}), 2, 700, true);
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(identity, replay.fragments, workspace);
    ufdsMaskConstructWithWorkspace(split, replay.fragments, workspace);
    assert(replay.fragments.size() == 4);
    for (size_t i = 1; i < replay.fragments.size(); i++)
        assert(replay.fragments[i].canonicalId == resolvedIds[i - 1]);
}

void testResidualCacheCanonicalIdWriteback()
{
    ufdsMaskWorkspace workspace(33, 32);
    const EdgeMask identity = makeMask({12, 13, 14, 15});
    const vector<int> identityIds =
        resolveFirstCacheReplay(identity, workspace);
    assert(identityIds.size() == 1);

    assemblyState identityHit;
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(
        identity,
        identityHit.fragments,
        workspace
    );
    assert(identityHit.fragments.size() == 1);
    assert(identityHit.fragments[0].canonicalId == identityIds[0]);
    assert(identityHit.fragments[0].connected);

    const EdgeMask split = makeMask({18, 19, 22, 23});
    const vector<int> splitIds = resolveFirstCacheReplay(split, workspace);
    assert(splitIds.size() == 2);

    assemblyState splitHit;
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(split, splitHit.fragments, workspace);
    assert(splitHit.fragments.size() == 2);
    assert(splitHit.fragments[0].canonicalId == splitIds[0]);
    assert(splitHit.fragments[1].canonicalId == splitIds[1]);
    assert(splitHit.fragments[0].connected);
    assert(splitHit.fragments[1].connected);
}

void testWideResidualCacheCanonicalIdWriteback()
{
    configurePathGraph(65);
    ufdsMaskWorkspace workspace(66, 65);
    const EdgeMask identity = makeMask({60, 61, 62, 63, 64});
    const vector<int> canonicalIds =
        resolveFirstCacheReplay(identity, workspace);
    assert(canonicalIds.size() == 1);

    assemblyState hit;
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(identity, hit.fragments, workspace);
    assert(hit.fragments.size() == 1);
    assert(hit.fragments[0].canonicalId == canonicalIds[0]);
    assert(hit.fragments[0].connected);

    const EdgeMask split = makeMask({0, 1, 63, 64});
    const vector<int> splitIds = resolveFirstCacheReplay(split, workspace);
    assert(splitIds.size() == 2);

    assemblyState splitHit;
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(split, splitHit.fragments, workspace);
    assert(splitHit.fragments.size() == 2);
    assert(splitHit.fragments[0].canonicalId == splitIds[0]);
    assert(splitHit.fragments[1].canonicalId == splitIds[1]);
}

int main()
{
    configurePathGraph(32);
    testMaskTrimmingMetadata();
    testFragmentMetadataAndSelectiveCanonisation();
    testResidualCacheCanonicalIdWriteback();
    testMultipleResidualCacheBindings();
    testWideResidualCacheCanonicalIdWriteback();

    bitsetHashTable.clear();
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(64);
    std::construct_at(std::addressof(allEdges));
}
