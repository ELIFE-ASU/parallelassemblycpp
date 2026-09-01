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
using IntegerVector = vector<int>;
using BooleanVector = vector<bool>;
using IntegerPair = pair<int, int>;

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
    assert(
        edgeCount < static_cast<size_t>(numeric_limits<int>::max()) &&
        edgeCount <= static_cast<size_t>(numeric_limits<unsigned int>::max())
    );
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    universeEdgeList.clear();
    targetMolecule = molGraph();
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(edgeCount);
    std::construct_at(std::addressof(allEdges));
    string atomType = "C";
    for (size_t atomIndex = 0; atomIndex <= edgeCount; atomIndex++)
        targetMolecule.addAtom(atomType);
    for (size_t edge = 0; edge < edgeCount; edge++)
    {
        const int firstVertex = static_cast<int>(edge);
        targetMolecule.addBond(firstVertex, firstVertex + 1, 1);
    }
    universeEdgeList = targetMolecule.writeEdgeList();
    prepareCanonicalisationGraph(targetMolecule, universeEdgeList);
    AtomMask::configure(edgeCount + 1);
    totalBonds = static_cast<unsigned int>(edgeCount);
}

/** Own one production duplicate level and keep its sealed spans alive. */
struct dagDuplicateLevelFixture
{
    dagDuplicateClassLevel level;
    duplicateClassIndexWorkspace classIndex;
    size_t duplicateSize;
    size_t fragmentCount;

    dagDuplicateLevelFixture(
        size_t duplicateSizeValue,
        size_t fragmentCountValue
    ):
        duplicateSize(duplicateSizeValue), fragmentCount(fragmentCountValue)
    {
        rebuild();
    }

    void rebuild()
    {
        level.reset(fragmentCount);
        classIndex.beginLevel();
    }

    void insert(int canonicalId, potentialDuplicate occurrence)
    {
        classIndex.getOrCreate(
            level,
            canonicalId,
            duplicateSize,
            fragmentCount
        ).insert(std::move(occurrence));
    }

    void seal()
    {
        level.seal();
    }

    duplicateClassEntry<dagDuplicateSet> &classEntry(int canonicalId)
    {
        const auto entry = find_if(
            level.classes.begin(),
            level.classes.end(),
            [&](const auto &candidate)
            {
                return candidate.canonicalId == canonicalId;
            }
        );
        assert(entry != level.classes.end());
        return *entry;
    }

    dagDuplicateSet &singleClass()
    {
        assert(level.classes.size() == 1);
        return level.classes.front().duplicates;
    }
};

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

    IntegerVector key;
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
    IntegerVector key;
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
    const IntegerVector firstKey = key;
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
    for (size_t index = 1; index < child.fragments.size(); index++)
        assert(child.fragments[index].canonicalId == unknownCanonicalId);

    IntegerVector key;
    assert(canoniseAssemblyStateAndBuildKey(child, key, workspace));
    vector<int> resolvedIds;
    for (size_t index = 1; index < child.fragments.size(); index++)
        resolvedIds.push_back(child.fragments[index].canonicalId);

    assemblyState replay;
    replay.appendFragment(makeMask({0, 1}), 2, 700, true);
    workspace.beginFragmentation();
    ufdsMaskConstructWithWorkspace(identity, replay.fragments, workspace);
    ufdsMaskConstructWithWorkspace(split, replay.fragments, workspace);
    assert(replay.fragments.size() == 4);
    for (size_t index = 1; index < replay.fragments.size(); index++)
        assert(
            replay.fragments[index].canonicalId == resolvedIds[index - 1]
        );
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

void testFragmentPairMatchingTraversal()
{
    dagDuplicateLevelFixture duplicatesStorage(2, 3);
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({0, 1}), 0, 0)
    );
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({2, 3}), 0, 2)
    );
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({1, 4}), 0, 5)
    );
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({8, 9}), 1, 1)
    );
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({10, 11}), 1, 4)
    );
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({16, 17}), 2, 3)
    );
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({18, 19}), 2, 6)
    );
    duplicatesStorage.insert(
        10,
        potentialDuplicate(makeMask({20, 21}), 2, 7)
    );
    duplicatesStorage.seal();
    dagDuplicateSet &duplicates = duplicatesStorage.singleClass();
    assert(duplicates.isValid());
    vector<int> fragmentsByOccurrence;
    for (const potentialDuplicate &occurrence : duplicates.list)
        fragmentsByOccurrence.push_back(occurrence.fragmentIndex);

    vector<pair<size_t, size_t>> expectedPairs;
    assert(duplicates.visitMatchingsInReverse(
        [&](validMatchings &matching, size_t first, size_t second)
        {
            assert(first < second);
            assert(
                matching.firstFragmentIndex == fragmentsByOccurrence[first]
            );
            assert(
                matching.secondFragmentIndex == fragmentsByOccurrence[second]
            );
            expectedPairs.emplace_back(
                matching.first.findFirst(),
                matching.second.findFirst()
            );
            return false;
        },
        [](validMatchings &) {return true;}
    ));
    assert(expectedPairs.size() == 27);
    sort(expectedPairs.begin(), expectedPairs.end());

    vector<pair<int, int>> visitedBlocks;
    vector<pair<size_t, size_t>> actualPairs;
    assert(!duplicates.hasDenseFragmentRuns());
    assert(duplicates.visitMatchingsByFragmentPairInReverse(
        [&](int firstFragment, int secondFragment, int selectedSize)
        {
            assert(firstFragment <= secondFragment);
            assert(static_cast<size_t>(selectedSize) == duplicates.size);
            visitedBlocks.emplace_back(
                firstFragment,
                secondFragment
            );
            return false;
        },
        [&](validMatchings &matching)
        {
            actualPairs.emplace_back(
                matching.first.findFirst(),
                matching.second.findFirst()
            );
            return true;
        }
    ));
    sort(actualPairs.begin(), actualPairs.end());
    assert(actualPairs == expectedPairs);
    sort(visitedBlocks.begin(), visitedBlocks.end());
    const vector<pair<int, int>> expectedBlocks = {
        {0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}
    };
    assert(visitedBlocks == expectedBlocks);

    size_t blockFilterCalls = 0;
    size_t matchingVisits = 0;
    assert(duplicates.visitMatchingsByFragmentPairInReverse(
        [&](int firstFragment, int secondFragment, int selectedSize)
        {
            ++blockFilterCalls;
            assert(static_cast<size_t>(selectedSize) == duplicates.size);
            return firstFragment == 0 && secondFragment == 2;
        },
        [&](validMatchings &)
        {
            ++matchingVisits;
            return true;
        }
    ));
    assert(blockFilterCalls == 6);
    assert(matchingVisits == expectedPairs.size() - 9);

    size_t stoppedVisits = 0;
    assert(!duplicates.visitMatchingsByFragmentPairInReverse(
        [](int, int, int) {return false;},
        [&](validMatchings &)
        {
            ++stoppedVisits;
            return stoppedVisits < 3;
        }
    ));
    assert(stoppedVisits == 3);

    // Empty diagonal blocks do not evaluate a fragment-pair filter.
    dagDuplicateLevelFixture singletonRunStorage(2, 2);
    for (size_t occurrence = 0; occurrence < 8; occurrence++)
    {
        singletonRunStorage.insert(
            20,
            potentialDuplicate(
                makeMask({2 * occurrence, 2 * occurrence + 1}),
                occurrence == 0 ? 0 : 1,
                static_cast<int>(occurrence)
            )
        );
    }
    singletonRunStorage.seal();
    dagDuplicateSet &singletonRun = singletonRunStorage.singleClass();
    assert(!singletonRun.hasDenseFragmentRuns());
    vector<pair<int, int>> singletonBlocks;
    size_t singletonVisits = 0;
    assert(singletonRun.visitMatchingsByFragmentPairInReverse(
        [&](int firstFragment, int secondFragment, int)
        {
            singletonBlocks.emplace_back(firstFragment, secondFragment);
            return false;
        },
        [&](validMatchings &)
        {
            ++singletonVisits;
            return true;
        }
    ));
    sort(singletonBlocks.begin(), singletonBlocks.end());
    const vector<pair<int, int>> expectedSingletonBlocks = {
        {0, 1}, {1, 1}
    };
    assert(singletonBlocks == expectedSingletonBlocks);
    assert(singletonVisits == 28);

    // Resume a diagonal after overlapping pairs precede its first match.
    dagDuplicateLevelFixture delayedDiagonalStorage(2, 1);
    delayedDiagonalStorage.insert(
        30,
        potentialDuplicate(makeMask({0, 1}), 0, 0)
    );
    delayedDiagonalStorage.insert(
        30,
        potentialDuplicate(makeMask({2, 3}), 0, 1)
    );
    delayedDiagonalStorage.insert(
        30,
        potentialDuplicate(makeMask({0, 4}), 0, 2)
    );
    delayedDiagonalStorage.insert(
        30,
        potentialDuplicate(makeMask({0, 5}), 0, 3)
    );
    delayedDiagonalStorage.seal();
    dagDuplicateSet &delayedDiagonal = delayedDiagonalStorage.singleClass();
    vector<pair<size_t, size_t>> expectedDelayedPairs;
    assert(delayedDiagonal.visitMatchingsInReverse(
        [&](validMatchings &matching, size_t, size_t)
        {
            expectedDelayedPairs.emplace_back(
                matching.first.findFirst(),
                matching.second.findFirst()
            );
            return false;
        },
        [](validMatchings &) {return true;}
    ));
    vector<pair<size_t, size_t>> actualDelayedPairs;
    size_t delayedFilterCalls = 0;
    assert(delayedDiagonal.visitMatchingsByFragmentPairInReverse(
        [&](int, int, int)
        {
            ++delayedFilterCalls;
            return false;
        },
        [&](validMatchings &matching)
        {
            actualDelayedPairs.emplace_back(
                matching.first.findFirst(),
                matching.second.findFirst()
            );
            return true;
        }
    ));
    sort(expectedDelayedPairs.begin(), expectedDelayedPairs.end());
    sort(actualDelayedPairs.begin(), actualDelayedPairs.end());
    assert(expectedDelayedPairs.size() == 3);
    assert(actualDelayedPairs == expectedDelayedPairs);
    assert(delayedFilterCalls == 1);

    dagDuplicateLevelFixture denseSmallStorage(2, 1);
    for (size_t occurrence = 0; occurrence < 8; occurrence++)
    {
        denseSmallStorage.insert(
            40,
            potentialDuplicate(
                makeMask({2 * occurrence, 2 * occurrence + 1}),
                0,
                static_cast<int>(occurrence)
            )
        );
    }
    denseSmallStorage.seal();
    dagDuplicateSet &denseSmall = denseSmallStorage.singleClass();
    assert(denseSmall.hasDenseFragmentRuns());

    // Sorted but sparse runs also retain the legacy traversal.
    dagDuplicateLevelFixture sparseDuplicatesStorage(2, 8);
    for (size_t occurrence = 0; occurrence < 16; occurrence++)
    {
        sparseDuplicatesStorage.insert(
            50,
            potentialDuplicate(
                makeMask({2 * occurrence, 2 * occurrence + 1}),
                static_cast<int>(occurrence / 2),
                static_cast<int>(occurrence)
            )
        );
    }
    sparseDuplicatesStorage.seal();
    dagDuplicateSet &sparseDuplicates =
        sparseDuplicatesStorage.singleClass();
    assert(!sparseDuplicates.hasDenseFragmentRuns());

    // Unexpected interleaving retains the legacy occurrence traversal. Keep
    // this case wide and sparse to cover that correctness fallback.
    const array<int, 16> interleavedFragments = {
        7, 0, 3, 7, 0, 3, 7, 3, 0, 7, 0, 3, 7, 0, 3, 7
    };
    dagDuplicateLevelFixture interleavedDuplicatesStorage(2, 8);
    for (size_t occurrence = 0;
         occurrence < interleavedFragments.size();
         occurrence++)
    {
        EdgeMask mask;
        mask.set(occurrence == 15 ? 0 : occurrence * 2);
        mask.set(occurrence * 2 + 1);
        interleavedDuplicatesStorage.insert(
            60,
            potentialDuplicate(
                std::move(mask),
                interleavedFragments[occurrence],
                static_cast<int>(occurrence)
            )
        );
    }
    interleavedDuplicatesStorage.seal();
    dagDuplicateSet &interleavedDuplicates =
        interleavedDuplicatesStorage.singleClass();
    assert(interleavedDuplicates.isValid());

    size_t expectedFallbackPairs = 0;
    assert(interleavedDuplicates.visitMatchingsInReverse(
        [&](validMatchings &matching, size_t first, size_t second)
        {
            assert(first < second);
            assert(
                matching.firstFragmentIndex == interleavedFragments[first]
            );
            assert(
                matching.secondFragmentIndex == interleavedFragments[second]
            );
            ++expectedFallbackPairs;
            return false;
        },
        [](validMatchings &) {return true;}
    ));
    assert(expectedFallbackPairs == 119);

    assert(!interleavedDuplicates.hasDenseFragmentRuns());
    size_t fallbackVisits = 0;
    assert(interleavedDuplicates.visitMatchingsInReverse(
        [&](validMatchings &matching, size_t first, size_t second)
        {
            assert(first < second);
            assert(
                matching.firstFragmentIndex == interleavedFragments[first]
            );
            assert(
                matching.secondFragmentIndex == interleavedFragments[second]
            );
            return false;
        },
        [&](validMatchings &)
        {
            ++fallbackVisits;
            return true;
        }
    ));
    assert(fallbackVisits == expectedFallbackPairs);
}

void testDuplicateClassCsrStorage()
{
    // Canonical-class packing is sorted by ID, but each class retains its
    // original insertion order when insertions are interleaved.
    dagDuplicateLevelFixture interleaved(2, 6);
    interleaved.insert(200, potentialDuplicate(makeMask({0, 1}), 0, 0));
    interleaved.insert(100, potentialDuplicate(makeMask({2, 3}), 1, 1));
    interleaved.insert(200, potentialDuplicate(makeMask({4, 5}), 2, 2));
    interleaved.insert(100, potentialDuplicate(makeMask({6, 7}), 3, 3));
    interleaved.insert(200, potentialDuplicate(makeMask({8, 9}), 4, 4));
    interleaved.insert(100, potentialDuplicate(makeMask({10, 11}), 5, 5));
    interleaved.seal();

    assert(interleaved.level.classes.size() == 2);
    assert(interleaved.level.classes[0].canonicalId == 100);
    assert(interleaved.level.classes[1].canonicalId == 200);
    const auto &firstClass = interleaved.level.classes[0].duplicates;
    const auto &secondClass = interleaved.level.classes[1].duplicates;
    assert(firstClass.list.size() == 3);
    assert(secondClass.list.size() == 3);
    assert(firstClass.list[0].duplicateIndex == 1);
    assert(firstClass.list[1].duplicateIndex == 3);
    assert(firstClass.list[2].duplicateIndex == 5);
    assert(secondClass.list[0].duplicateIndex == 0);
    assert(secondClass.list[1].duplicateIndex == 2);
    assert(secondClass.list[2].duplicateIndex == 4);

    // Sparse rows retain the complete fragment domain while storing only
    // non-empty fragment unions.
    dagDuplicateLevelFixture sparse(2, 12);
    sparse.insert(300, potentialDuplicate(makeMask({0, 1}), 1, 0));
    sparse.insert(300, potentialDuplicate(makeMask({64, 65}), 9, 1));
    sparse.insert(300, potentialDuplicate(makeMask({65, 66}), 9, 2));
    sparse.seal();
    auto &sparseEntry = sparse.classEntry(300);
    const duplicateFragmentMaskList sparseMasks =
        sparse.level.fragmentMasks(sparseEntry);
    assert(sparse.level.fragmentCount() == 12);
    assert(sparseEntry.duplicates.fragmentCount == 12);
    assert(sparseMasks.fragmentCount == 12);
    assert(sparseMasks.size() == 2);
    assert(sparseMasks[0].fragment == 1);
    assert(sparseMasks[1].fragment == 9);
    assert(sparseMasks.maskCount(0) == 0);
    assert(sparseMasks.maskCount(1) == 2);
    assert(sparseMasks.maskCount(8) == 0);
    assert(sparseMasks.maskCount(9) == 3);
    assert(sparseMasks.maskCount(11) == 0);

    // The same owned storage can be rebuilt after both sealed and unsealed
    // populations without retaining old classes, occurrences, or mask rows.
    dagDuplicateLevelFixture reusable(2, 5);
    reusable.insert(400, potentialDuplicate(makeMask({12, 13}), 0, 0));
    reusable.insert(400, potentialDuplicate(makeMask({14, 15}), 1, 1));
    reusable.seal();
    assert(reusable.level.classes.size() == 1);
    assert(reusable.singleClass().list.size() == 2);

    reusable.rebuild();
    assert(reusable.level.empty());
    assert(reusable.level.fragmentCount() == 5);
    reusable.insert(500, potentialDuplicate(makeMask({16, 17}), 2, 2));
    reusable.insert(500, potentialDuplicate(makeMask({18, 19}), 3, 3));

    // Reset the populated but unsealed staging representation as well.
    reusable.rebuild();
    assert(reusable.level.empty());
    assert(reusable.level.fragmentCount() == 5);
    reusable.insert(600, potentialDuplicate(makeMask({20, 21}), 4, 4));
    reusable.insert(600, potentialDuplicate(makeMask({22, 23}), 4, 5));
    reusable.seal();
    assert(reusable.level.classes.size() == 1);
    assert(reusable.level.classes.front().canonicalId == 600);
    assert(reusable.singleClass().list.size() == 2);
    const duplicateFragmentMaskList rebuiltMasks =
        reusable.level.fragmentMasks(reusable.level.classes.front());
    assert(rebuiltMasks.fragmentCount == 5);
    assert(rebuiltMasks.size() == 1);
    assert(rebuiltMasks[0].fragment == 4);
    assert(rebuiltMasks.maskCount(4) == 4);
}

int main()
{
    // Exercise grouped runs and their fallback with a multiword mask domain.
    configurePathGraph(96);
    testDuplicateClassCsrStorage();
    testFragmentPairMatchingTraversal();
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
