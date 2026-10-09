// This executable is also built by Release/CI presets; keep its checks active.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <csignal>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;
using IntegerVector = vector<int>;
using BooleanVector = vector<bool>;
using IntegerPair = pair<int, int>;

#include "../src/activeWordMask.h"

constexpr int ceilLog2(int value)
{
    return std::bit_width(static_cast<unsigned int>(value - 1));
}

#include "../src/globalPrimitives.h"
#include "../src/ufds.h"
#include "../src/molGraph.h"
#include "../src/treeCanon.h"

#include "graphTestFixtures.h"

using graphTestFixtures::edgeSpec;
using graphTestFixtures::makeGraph;

treeCanonForm canonicalForm(molGraph &tree)
{
    return centroidTreeCanon(tree, 0);
}

bool equivalentTrees(molGraph left, molGraph right)
{
    clearTreeCanonInterner();
    const treeCanonForm leftForm = canonicalForm(left);
    const treeCanonForm rightForm = canonicalForm(right);
    return leftForm == rightForm;
}

vector<int> reversePermutation(size_t vertexCount)
{
    assert(vertexCount <= static_cast<size_t>(numeric_limits<int>::max()));
    vector<int> permutation(vertexCount);
    for (size_t index = 0; index < vertexCount; index++)
        permutation[index] = static_cast<int>(vertexCount - index - 1);
    return permutation;
}

void testCentroidForms()
{
    molGraph singleton = makeGraph({"C"}, {});
    clearTreeCanonInterner();
    const treeCanonForm singletonForm = canonicalForm(singleton);
    assert(!singletonForm.empty());
    assert(singletonForm.second == 0);

    molGraph pair = makeGraph({"C", "N"}, {{0, 1, 2}});
    clearTreeCanonInterner();
    const treeCanonForm pairForm = canonicalForm(pair);
    assert(!pairForm.empty());
    assert(pairForm.second != 0);
    assert(pairForm.centralBond == 2);
}

void testLabelsAndPermutation()
{
    const vector<string> labels{"C", "N", "O", "C", "S", "C", "F", "C"};
    const vector<edgeSpec> edges{
        {0, 1, 1}, {1, 2, 2}, {1, 3, 1}, {3, 4, 1},
        {4, 5, 2}, {4, 6, 1}, {6, 7, 1}
    };
    const vector<int> permutation{5, 1, 7, 0, 6, 3, 2, 4};
    assert(equivalentTrees(
        makeGraph(labels, edges),
        makeGraph(labels, edges, permutation, true)
    ));

    vector<string> movedAtom = labels;
    swap(movedAtom[2], movedAtom[3]);
    assert(!equivalentTrees(
        makeGraph(labels, edges),
        makeGraph(movedAtom, edges)
    ));
}

void testBondPlacement()
{
    const vector<string> labels(5, "C");
    const vector<edgeSpec> terminalDouble{
        {0, 1, 2}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}
    };
    const vector<edgeSpec> centralDouble{
        {0, 1, 1}, {1, 2, 2}, {2, 3, 1}, {3, 4, 1}
    };
    assert(!equivalentTrees(
        makeGraph(labels, terminalDouble),
        makeGraph(labels, centralDouble)
    ));
}

void testWideBondLabels()
{
    molGraph lowLabel = makeGraph({"C", "N"}, {{0, 1, 1}});
    molGraph highLabel = makeGraph({"C", "N"}, {{0, 1, 257}});
    assert(!equivalentTrees(lowLabel, highLabel));

    clearTreeCanonInterner();
    assert(canonicalForm(highLabel).centralBond == 257);

    molGraph uniqueEdges = makeGraph(
        {"C", "N", "C", "N"},
        {{0, 1, 1}, {2, 3, 257}}
    );
    vector<MoleculeEdge> removedEdges;
    const molGraph processed = preprocessWriteback(uniqueEdges, removedEdges);
    assert(removedEdges.size() == 2);
    assert(processed.totalBonds == 0);
}

void testSameDegreeNonIsomorphs()
{
    const vector<string> labels(7, "C");
    const vector<edgeSpec> adjacentBranches{
        {0, 1, 1}, {0, 2, 1}, {2, 3, 1},
        {0, 4, 1}, {1, 5, 1}, {1, 6, 1}
    };
    const vector<edgeSpec> separatedBranches{
        {0, 2, 1}, {2, 1, 1}, {0, 3, 1},
        {0, 4, 1}, {1, 5, 1}, {1, 6, 1}
    };
    assert(!equivalentTrees(
        makeGraph(labels, adjacentBranches),
        makeGraph(labels, separatedBranches)
    ));
}

void testLongPaths()
{
    for (const size_t nodeCount : {511U, 512U, 32766U, 32767U})
    {
        vector<string> labels(nodeCount, "C");
        vector<edgeSpec> edges;
        edges.reserve(nodeCount - 1);
        for (size_t node = 1; node < nodeCount; node++)
        {
            labels[node] = node % 7 == 0 ? "N" : "C";
            edges.emplace_back(
                static_cast<int>(node - 1),
                static_cast<int>(node),
                static_cast<short>(node % 5 == 0 ? 2 : 1)
            );
        }
        assert(equivalentTrees(
            makeGraph(labels, edges),
            makeGraph(labels, edges, reversePermutation(nodeCount), true)
        ));
    }
}

void testHighDegreeTree()
{
    constexpr size_t nodeCount = 512;
    vector<string> labels(nodeCount, "C");
    labels[0] = "N";
    vector<edgeSpec> edges;
    edges.reserve(nodeCount - 1);
    for (size_t node = 1; node < nodeCount; node++)
    {
        labels[node] = node % 3 == 0 ? "O" : "C";
        edges.emplace_back(
            0,
            static_cast<int>(node),
            static_cast<short>(node % 11 == 0 ? 2 : 1)
        );
    }
    vector<int> permutation(nodeCount);
    for (size_t node = 0; node < nodeCount; node++)
        permutation[node] = static_cast<int>((node * 173) % nodeCount);
    assert(equivalentTrees(
        makeGraph(labels, edges),
        makeGraph(labels, edges, permutation, true)
    ));
}

void testUnsupportedGraphsFallBack()
{
    molGraph cycle = makeGraph(
        {"C", "C", "C"},
        {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}}
    );
    clearTreeCanonInterner();
    assert(canonicalForm(cycle).empty());

    molGraph forest = makeGraph(
        {"C", "C", "C", "C"},
        {{0, 1, 1}, {1, 2, 1}}
    );
    clearTreeCanonInterner();
    assert(canonicalForm(forest).empty());

    molGraph sentinel = makeGraph(
        {"C", "X", "C"},
        {{0, 1, 1}, {1, 2, 1}}
    );
    clearTreeCanonInterner();
    assert(canonicalForm(sentinel).empty());
}

void testWideUfdsSplit()
{
    EdgeMask::configure(513);
    {
        ufdsSplit splitter;
        splitter.elements.resize(514);
        splitter.reset();

        // Insert the high component first and include atom/edge bit 512. The
        // splitter must still emit components in atom order and reconstruct
        // the complete nine-word masks without the former fixed arrays.
        splitter.doubleInsert(512, 513, 512);
        splitter.insert(511, 512, 511);
        splitter.doubleInsert(0, 1, 0);
        splitter.insert(2, 1, 1);

        vector<assemblyFragment> fragments;
        splitter.splitWithBuffers(fragments);

        assert(fragments.size() == 2);
        assert(fragments[0].edgeCount == 2);
        assert(fragments[1].edgeCount == 2);
        assert(fragments[0].connected);
        assert(fragments[1].connected);
        assert(fragments[0].mask.count() == 2);
        assert(fragments[0].mask.test(0));
        assert(fragments[0].mask.test(1));
        assert(fragments[1].mask.count() == 2);
        assert(fragments[1].mask.test(511));
        assert(fragments[1].mask.test(512));

        // Reuse the workspace to verify reset clears the sparse wide-word
        // tracking rather than leaking the preceding high component.
        splitter.reset();
        fragments.clear();
        splitter.doubleInsert(512, 513, 7);
        splitter.insert(511, 512, 8);
        splitter.splitWithBuffers(fragments);
        assert(fragments.size() == 1);
        assert(fragments[0].edgeCount == 2);
        assert(fragments[0].connected);
        assert(fragments[0].mask.count() == 2);
        assert(fragments[0].mask.test(7));
        assert(fragments[0].mask.test(8));
        assert(!fragments[0].mask.test(511));
        assert(!fragments[0].mask.test(512));
    }
    EdgeMask::configure(64);
}

int main()
{
    testCentroidForms();
    testLabelsAndPermutation();
    testBondPlacement();
    testWideBondLabels();
    testSameDegreeNonIsomorphs();
    testLongPaths();
    testHighDegreeTree();
    testUnsupportedGraphsFallBack();
    testWideUfdsSplit();
    clearTreeCanonInterner();
    assert(treeCanonInterner.empty());
    assert(treeCanonAtomInterner.empty());
    return 0;
}
