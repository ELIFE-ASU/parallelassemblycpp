#include <algorithm>
#include <array>
#include <atomic>
#include <barrier>
#include <bit>
#include <cassert>
#include <cstdint>
#include <csignal>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;
using vi = vector<int>;
using vb = vector<bool>;
using pii = pair<int, int>;

// Model the worker-local state used by the OpenMP executable. Each test thread
// must have its own mask arena, canonical caches, interners, and scratch.
#define ASSEMBLYCPP_SEARCH_LOCAL thread_local
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
#include "../v5/graphHashes.h"

using edgeSpec = std::tuple<int, int, short>;

struct namedGraph
{
    string name;
    molGraph graph;
};

[[noreturn]] void fail(const string &message)
{
    cerr << "cyclic canonicalisation test failed: " << message << '\n';
    abort();
}

void require(bool condition, const string &message)
{
    if (!condition) fail(message);
}

vector<int> shuffledPermutation(size_t size, uint32_t seed)
{
    vector<int> permutation(size);
    iota(permutation.begin(), permutation.end(), 0);
    for (size_t remaining = size; remaining > 1; remaining--)
    {
        seed = seed * 1664525U + 1013904223U;
        swap(permutation[remaining - 1], permutation[seed % remaining]);
    }
    return permutation;
}

molGraph makeGraph(
    const vector<string> &labels,
    const vector<edgeSpec> &edges,
    const vector<int> &oldToNew = {},
    bool reverseEdges = false
)
{
    vector<int> permutation = oldToNew;
    if (permutation.empty())
    {
        permutation.resize(labels.size());
        iota(permutation.begin(), permutation.end(), 0);
    }
    require(permutation.size() == labels.size(), "invalid permutation size");

    vector<unsigned char> seen(labels.size(), 0);
    vector<string> permutedLabels(labels.size());
    for (size_t old = 0; old < labels.size(); old++)
    {
        const int replacement = permutation[old];
        require(
            replacement >= 0 &&
                static_cast<size_t>(replacement) < labels.size() &&
                !seen[replacement],
            "permutation is not a bijection"
        );
        seen[replacement] = 1;
        permutedLabels[replacement] = labels[old];
    }

    molGraph result;
    for (string &label : permutedLabels) result.addAtom(label);
    const auto addOne = [&](const edgeSpec &edge)
    {
        const auto [left, right, bondType] = edge;
        result.addBond(permutation[left], permutation[right], bondType);
    };
    if (reverseEdges)
    {
        for (auto edge = edges.rbegin(); edge != edges.rend(); ++edge)
            addOne(*edge);
    }
    else
    {
        for (const edgeSpec &edge : edges) addOne(edge);
    }
    return result;
}

/**
 * Two copies each of K3,3 and the triangular prism. Both graph families have
 * six equally-labelled vertices, nine equal bonds, and degree three at every
 * vertex, so the cheap shared invariant deliberately collides. They are not
 * isomorphic, while copies within one family are exact equivalents.
 */
molGraph makeSharedRegistryCollisionGraph()
{
    vector<string> labels(24, "C");
    vector<edgeSpec> edges;
    edges.reserve(36);
    const auto addCompleteBipartite = [&](int base)
    {
        for (int left = 0; left < 3; ++left)
        {
            for (int right = 3; right < 6; ++right)
                edges.emplace_back(base + left, base + right, 1);
        }
    };
    const auto addTriangularPrism = [&](int base)
    {
        edges.emplace_back(base, base + 1, 1);
        edges.emplace_back(base + 1, base + 2, 1);
        edges.emplace_back(base + 2, base, 1);
        edges.emplace_back(base + 3, base + 4, 1);
        edges.emplace_back(base + 4, base + 5, 1);
        edges.emplace_back(base + 5, base + 3, 1);
        edges.emplace_back(base, base + 3, 1);
        edges.emplace_back(base + 1, base + 4, 1);
        edges.emplace_back(base + 2, base + 5, 1);
    };

    addCompleteBipartite(0);
    addCompleteBipartite(6);
    addTriangularPrism(12);
    addTriangularPrism(18);
    return makeGraph(labels, edges);
}

/**
 * Two copies each of non-isomorphic trees with the same degree sequence,
 * followed by two triangles carrying the same three-edge pendant chain, one
 * labelled seed edge, and one cyclic producer-only seed class.
 */
molGraph makeSharedRegistryTreeInternerGraph()
{
    vector<string> labels(41, "C");
    labels[36] = "N";
    labels[37] = "O";
    labels[38] = "P";
    labels[39] = "P";
    labels[40] = "P";
    vector<edgeSpec> edges;
    edges.reserve(33);
    const auto addForkedTree = [&](int base)
    {
        edges.emplace_back(base, base + 1, 1);
        edges.emplace_back(base, base + 2, 1);
        edges.emplace_back(base, base + 3, 1);
        edges.emplace_back(base + 1, base + 4, 1);
        edges.emplace_back(base + 2, base + 5, 1);
    };
    const auto addBentArmTree = [&](int base)
    {
        edges.emplace_back(base, base + 1, 1);
        edges.emplace_back(base, base + 2, 1);
        edges.emplace_back(base, base + 3, 1);
        edges.emplace_back(base + 1, base + 4, 1);
        edges.emplace_back(base + 4, base + 5, 1);
    };
    const auto addCycleWithPendantChain = [&](int base)
    {
        edges.emplace_back(base, base + 1, 1);
        edges.emplace_back(base + 1, base + 2, 1);
        edges.emplace_back(base + 2, base, 1);
        edges.emplace_back(base, base + 3, 1);
        edges.emplace_back(base + 3, base + 4, 1);
        edges.emplace_back(base + 4, base + 5, 1);
    };

    addForkedTree(0);
    addForkedTree(6);
    addBentArmTree(12);
    addBentArmTree(18);
    addCycleWithPendantChain(24);
    addCycleWithPendantChain(30);
    edges.emplace_back(36, 37, 2);
    edges.emplace_back(38, 39, 3);
    edges.emplace_back(39, 40, 3);
    edges.emplace_back(40, 38, 3);
    return makeGraph(labels, edges);
}

using bondLabelBag = vector<short>;
using neighbourInvariant = tuple<string, bondLabelBag, bondLabelBag>;

struct exactGraph
{
    vector<string> labels;
    vector<vector<bondLabelBag>> bonds;
    vector<bondLabelBag> incidentBondLabels;
    vector<vector<neighbourInvariant>> neighbourhoods;
};

exactGraph describeGraph(const molGraph &graph)
{
    const size_t size = graph.mg.size();
    exactGraph result;
    result.labels.reserve(size);
    result.bonds.assign(size, vector<bondLabelBag>(size));
    result.incidentBondLabels.resize(size);
    result.neighbourhoods.resize(size);

    for (const atom &vertex : graph.mg) result.labels.push_back(vertex.type);
    for (size_t first = 0; first < size; first++)
    {
        for (const bond &edge : graph.mg[first].list)
        {
            require(
                edge.n >= 0 && static_cast<size_t>(edge.n) < size,
                "graph contains an out-of-range bond endpoint"
            );
            result.incidentBondLabels[first].push_back(edge.type);

            const size_t second = static_cast<size_t>(edge.n);
            if (first < second)
            {
                result.bonds[first][second].push_back(edge.type);
                result.bonds[second][first].push_back(edge.type);
            }
            else if (first == second)
            {
                // molGraph records both ends of a self-loop in the same list.
                // Keeping both entries still gives an exact, consistent
                // multiplicity for comparisons between molGraph instances.
                result.bonds[first][first].push_back(edge.type);
            }
        }
        sort(
            result.incidentBondLabels[first].begin(),
            result.incidentBondLabels[first].end()
        );
    }

    for (size_t first = 0; first < size; first++)
    {
        for (size_t second = 0; second < size; second++)
        {
            bondLabelBag &labels = result.bonds[first][second];
            sort(labels.begin(), labels.end());
            if (labels.empty()) continue;
            result.neighbourhoods[first].emplace_back(
                result.labels[second],
                result.incidentBondLabels[second],
                labels
            );
        }
        sort(
            result.neighbourhoods[first].begin(),
            result.neighbourhoods[first].end()
        );
    }
    return result;
}

bool sameVertexInvariant(
    const exactGraph &left,
    size_t leftVertex,
    const exactGraph &right,
    size_t rightVertex
)
{
    return
        left.labels[leftVertex] == right.labels[rightVertex] &&
        left.incidentBondLabels[leftVertex] ==
            right.incidentBondLabels[rightVertex] &&
        left.bonds[leftVertex][leftVertex] ==
            right.bonds[rightVertex][rightVertex] &&
        left.neighbourhoods[leftVertex] == right.neighbourhoods[rightVertex];
}

struct exactMatcher
{
    const exactGraph &left;
    const exactGraph &right;
    vector<vector<size_t>> candidates;
    vector<int> leftToRight;
    vector<unsigned char> rightUsed;

    bool pairCompatible(size_t leftVertex, size_t rightVertex) const
    {
        for (size_t mappedLeft = 0; mappedLeft < leftToRight.size(); mappedLeft++)
        {
            if (leftToRight[mappedLeft] < 0) continue;
            const size_t mappedRight = static_cast<size_t>(leftToRight[mappedLeft]);
            if (
                left.bonds[leftVertex][mappedLeft] !=
                right.bonds[rightVertex][mappedRight]
            )
            {
                return false;
            }
        }
        return true;
    }

    bool search(size_t matched)
    {
        if (matched == leftToRight.size()) return true;

        size_t bestLeft = leftToRight.size();
        vector<size_t> bestRights;
        size_t bestCount = numeric_limits<size_t>::max();

        for (size_t leftVertex = 0; leftVertex < leftToRight.size(); leftVertex++)
        {
            if (leftToRight[leftVertex] >= 0) continue;

            vector<size_t> available;
            for (size_t rightVertex : candidates[leftVertex])
            {
                if (
                    !rightUsed[rightVertex] &&
                    pairCompatible(leftVertex, rightVertex)
                )
                {
                    available.push_back(rightVertex);
                }
            }
            if (available.empty()) return false;
            if (available.size() < bestCount)
            {
                bestLeft = leftVertex;
                bestRights = std::move(available);
                bestCount = bestRights.size();
            }
        }

        require(bestLeft < leftToRight.size(), "exact matcher lost an unmapped vertex");
        for (size_t rightVertex : bestRights)
        {
            leftToRight[bestLeft] = static_cast<int>(rightVertex);
            rightUsed[rightVertex] = 1;
            if (search(matched + 1)) return true;
            rightUsed[rightVertex] = 0;
            leftToRight[bestLeft] = -1;
        }
        return false;
    }
};

bool exactlyIsomorphic(const molGraph &leftInput, const molGraph &rightInput)
{
    if (leftInput.mg.size() != rightInput.mg.size()) return false;

    const exactGraph left = describeGraph(leftInput);
    const exactGraph right = describeGraph(rightInput);
    const size_t size = left.labels.size();
    vector<vector<size_t>> candidates(size);
    for (size_t leftVertex = 0; leftVertex < size; leftVertex++)
    {
        for (size_t rightVertex = 0; rightVertex < size; rightVertex++)
        {
            if (sameVertexInvariant(left, leftVertex, right, rightVertex))
                candidates[leftVertex].push_back(rightVertex);
        }
        if (candidates[leftVertex].empty()) return false;
    }

    exactMatcher matcher{
        left,
        right,
        std::move(candidates),
        vector<int>(size, -1),
        vector<unsigned char>(size, 0)
    };
    return matcher.search(0);
}

vector<edgeSpec> cycleEdges(int size, short bondType = 1)
{
    vector<edgeSpec> edges;
    edges.reserve(size);
    for (int vertex = 0; vertex < size; vertex++)
        edges.emplace_back(vertex, (vertex + 1) % size, bondType);
    return edges;
}

vector<edgeSpec> triangularPrismEdges()
{
    return {
        {0, 1, 1}, {1, 2, 1}, {2, 0, 1},
        {3, 4, 1}, {4, 5, 1}, {5, 3, 1},
        {0, 3, 1}, {1, 4, 1}, {2, 5, 1}
    };
}

vector<edgeSpec> completeBipartite33Edges()
{
    vector<edgeSpec> edges;
    for (int left = 0; left < 3; left++)
        for (int right = 3; right < 6; right++)
            edges.emplace_back(left, right, 1);
    return edges;
}

vector<edgeSpec> pentagonalPrismEdges()
{
    vector<edgeSpec> edges;
    for (int vertex = 0; vertex < 5; vertex++)
    {
        edges.emplace_back(vertex, (vertex + 1) % 5, 1);
        edges.emplace_back(vertex + 5, (vertex + 1) % 5 + 5, 1);
        edges.emplace_back(vertex, vertex + 5, 1);
    }
    return edges;
}

vector<edgeSpec> petersenEdges()
{
    vector<edgeSpec> edges;
    for (int vertex = 0; vertex < 5; vertex++)
    {
        edges.emplace_back(vertex, (vertex + 1) % 5, 1);
        edges.emplace_back(vertex, vertex + 5, 1);
    }
    edges.insert(
        edges.end(),
        {
            {5, 7, 1}, {7, 9, 1}, {9, 6, 1},
            {6, 8, 1}, {8, 5, 1}
        }
    );
    return edges;
}

vector<edgeSpec> rookGraphEdges()
{
    vector<edgeSpec> edges;
    for (int first = 0; first < 16; first++)
    {
        for (int second = first + 1; second < 16; second++)
        {
            if (first / 4 == second / 4 || first % 4 == second % 4)
                edges.emplace_back(first, second, 1);
        }
    }
    return edges;
}

vector<edgeSpec> shrikhandeGraphEdges()
{
    constexpr array<pair<int, int>, 6> offsets{{
        {1, 0}, {3, 0}, {0, 1}, {0, 3}, {1, 1}, {3, 3}
    }};
    set<pair<int, int>> uniqueEdges;
    for (int row = 0; row < 4; row++)
    {
        for (int column = 0; column < 4; column++)
        {
            const int first = row * 4 + column;
            for (const auto &[rowOffset, columnOffset] : offsets)
            {
                const int second = ((row + rowOffset) % 4) * 4 +
                    (column + columnOffset) % 4;
                uniqueEdges.emplace(min(first, second), max(first, second));
            }
        }
    }
    vector<edgeSpec> edges;
    for (const auto &[first, second] : uniqueEdges)
        edges.emplace_back(first, second, 1);
    return edges;
}

void appendPermutedPair(
    vector<namedGraph> &corpus,
    const string &name,
    const vector<string> &labels,
    const vector<edgeSpec> &edges,
    uint32_t seed
)
{
    corpus.push_back({name, makeGraph(labels, edges)});
    corpus.push_back({
        name + "-permuted",
        makeGraph(labels, edges, shuffledPermutation(labels.size(), seed), true)
    });
}

vector<namedGraph> makeDifferentialCorpus()
{
    vector<namedGraph> corpus;
    appendPermutedPair(corpus, "triangle", vector<string>(3, "C"), cycleEdges(3), 3);
    appendPermutedPair(corpus, "square", vector<string>(4, "C"), cycleEdges(4), 5);
    appendPermutedPair(corpus, "pentagon", vector<string>(5, "C"), cycleEdges(5), 7);

    vector<string> adjacentNitrogens(6, "C");
    adjacentNitrogens[0] = adjacentNitrogens[1] = "N";
    vector<string> oppositeNitrogens(6, "C");
    oppositeNitrogens[0] = oppositeNitrogens[3] = "N";
    appendPermutedPair(
        corpus, "adjacent-nitrogens", adjacentNitrogens, cycleEdges(6), 11
    );
    appendPermutedPair(
        corpus, "opposite-nitrogens", oppositeNitrogens, cycleEdges(6), 13
    );

    vector<edgeSpec> adjacentDoubleBonds = cycleEdges(6);
    std::get<2>(adjacentDoubleBonds[0]) = 2;
    std::get<2>(adjacentDoubleBonds[1]) = 2;
    vector<edgeSpec> separatedDoubleBonds = cycleEdges(6);
    std::get<2>(separatedDoubleBonds[0]) = 2;
    std::get<2>(separatedDoubleBonds[3]) = 2;
    appendPermutedPair(
        corpus,
        "adjacent-double-bonds",
        vector<string>(6, "C"),
        adjacentDoubleBonds,
        17
    );
    appendPermutedPair(
        corpus,
        "opposite-double-bonds",
        vector<string>(6, "C"),
        separatedDoubleBonds,
        19
    );

    vector<edgeSpec> adjacentAttachments = cycleEdges(4);
    adjacentAttachments.insert(
        adjacentAttachments.end(),
        {{0, 4, 1}, {4, 5, 2}, {1, 6, 1}, {6, 7, 2}}
    );
    vector<edgeSpec> oppositeAttachments = cycleEdges(4);
    oppositeAttachments.insert(
        oppositeAttachments.end(),
        {{0, 4, 1}, {4, 5, 2}, {2, 6, 1}, {6, 7, 2}}
    );
    vector<string> attachmentLabels(8, "C");
    attachmentLabels[5] = attachmentLabels[7] = "O";
    appendPermutedPair(
        corpus,
        "adjacent-attached-trees",
        attachmentLabels,
        adjacentAttachments,
        23
    );
    appendPermutedPair(
        corpus,
        "opposite-attached-trees",
        attachmentLabels,
        oppositeAttachments,
        29
    );

    appendPermutedPair(
        corpus,
        "triangular-prism",
        vector<string>(6, "C"),
        triangularPrismEdges(),
        31
    );
    appendPermutedPair(
        corpus,
        "complete-bipartite-3-3",
        vector<string>(6, "C"),
        completeBipartite33Edges(),
        37
    );
    appendPermutedPair(
        corpus,
        "pentagonal-prism",
        vector<string>(10, "C"),
        pentagonalPrismEdges(),
        41
    );
    appendPermutedPair(
        corpus,
        "petersen",
        vector<string>(10, "C"),
        petersenEdges(),
        43
    );
    appendPermutedPair(
        corpus,
        "rook-4x4",
        vector<string>(16, "C"),
        rookGraphEdges(),
        47
    );
    appendPermutedPair(
        corpus,
        "shrikhande",
        vector<string>(16, "C"),
        shrikhandeGraphEdges(),
        53
    );

    vector<edgeSpec> twoTriangles = cycleEdges(3);
    vector<edgeSpec> secondTriangle{{3, 4, 1}, {4, 5, 1}, {5, 3, 1}};
    twoTriangles.insert(twoTriangles.end(), secondTriangle.begin(), secondTriangle.end());
    appendPermutedPair(
        corpus,
        "two-disconnected-triangles",
        vector<string>(6, "C"),
        twoTriangles,
        59
    );
    appendPermutedPair(
        corpus,
        "six-cycle",
        vector<string>(6, "C"),
        cycleEdges(6),
        61
    );

    vector<edgeSpec> cyclicAndTree = cycleEdges(3);
    cyclicAndTree.insert(cyclicAndTree.end(), {{3, 4, 1}, {4, 5, 2}});
    vector<string> cyclicAndTreeLabels(6, "C");
    cyclicAndTreeLabels[5] = "N";
    appendPermutedPair(
        corpus,
        "disconnected-cycle-and-tree",
        cyclicAndTreeLabels,
        cyclicAndTree,
        67
    );

    appendPermutedPair(
        corpus,
        "forest-fallback",
        vector<string>{"C", "N", "C", "N", "O"},
        vector<edgeSpec>{{0, 1, 1}, {1, 2, 2}, {3, 4, 1}},
        71
    );
    appendPermutedPair(corpus, "empty-fallback", {}, {}, 73);
    appendPermutedPair(corpus, "singleton-fallback", {"C"}, {}, 79);

    appendPermutedPair(
        corpus,
        "labelled-parallel-edges",
        vector<string>{"C", "N"},
        vector<edgeSpec>{{0, 1, 1}, {0, 1, 2}},
        83
    );
    appendPermutedPair(
        corpus,
        "wide-bond-label-low",
        vector<string>(3, "C"),
        vector<edgeSpec>{{0, 1, 1}, {1, 2, 1}, {2, 0, 1}},
        89
    );
    appendPermutedPair(
        corpus,
        "wide-bond-label-high",
        vector<string>(3, "C"),
        vector<edgeSpec>{{0, 1, 257}, {1, 2, 1}, {2, 0, 1}},
        97
    );
    return corpus;
}

void testAgainstExactMatcher()
{
    clearTreeCanonInterner();
    vector<namedGraph> corpus = makeDifferentialCorpus();
    vector<cyclicCanonForm> forms;
    forms.reserve(corpus.size());
    for (namedGraph &entry : corpus)
    {
        forms.push_back(canonicaliseCyclicGraph(entry.graph));
        require(
            entry.graph.mg.empty() || !forms.back().empty(),
            entry.name + " produced an empty form"
        );
    }

    for (size_t left = 0; left < corpus.size(); left++)
    {
        for (size_t right = left; right < corpus.size(); right++)
        {
            const bool canonicalEquivalent = forms[left] == forms[right];
            const bool exactlyEquivalent = exactlyIsomorphic(
                corpus[left].graph,
                corpus[right].graph
            );
            if (canonicalEquivalent != exactlyEquivalent)
            {
                fail(
                    corpus[left].name + " vs " + corpus[right].name +
                    ": canonical=" + to_string(canonicalEquivalent) +
                    ", exact=" + to_string(exactlyEquivalent)
                );
            }
            if (canonicalEquivalent)
            {
                require(
                    forms[left].hash() == forms[right].hash(),
                    corpus[left].name + " vs " + corpus[right].name +
                        " have equal codes but unequal hashes"
                );
            }
        }
    }
}

void testLargeAttachedTrees()
{
    clearTreeCanonInterner();
    vector<string> labels(132, "C");
    labels.back() = "N";
    vector<edgeSpec> edges = cycleEdges(3);
    for (int vertex = 3; vertex < 132; vertex++)
        edges.emplace_back(vertex == 3 ? 0 : vertex - 1, vertex, vertex % 11 == 0 ? 2 : 1);

    molGraph original = makeGraph(labels, edges);
    molGraph permuted = makeGraph(
        labels,
        edges,
        shuffledPermutation(labels.size(), 89),
        true
    );
    const cyclicCanonForm originalForm = canonicaliseCyclicGraph(original);
    const cyclicCanonForm permutedForm = canonicaliseCyclicGraph(permuted);
    require(originalForm == permutedForm, "long attached tree is permutation-dependent");
    require(
        originalForm.hash() == permutedForm.hash(),
        "long attached tree hashes differ"
    );
    require(
        exactlyIsomorphic(original, permuted),
        "exact matcher rejected long-tree permutation"
    );

    vector<string> highDegreeLabels(68, "C");
    vector<edgeSpec> highDegreeEdges = cycleEdges(3);
    for (int vertex = 3; vertex < 68; vertex++)
        highDegreeEdges.emplace_back(0, vertex, vertex % 7 == 0 ? 2 : 1);
    molGraph highDegree = makeGraph(highDegreeLabels, highDegreeEdges);
    molGraph highDegreePermuted = makeGraph(
        highDegreeLabels,
        highDegreeEdges,
        shuffledPermutation(highDegreeLabels.size(), 97),
        true
    );
    require(
        canonicaliseCyclicGraph(highDegree) ==
            canonicaliseCyclicGraph(highDegreePermuted),
        "high-degree attached leaves are permutation-dependent"
    );
}

void testCachedGraphHashIsSelfContained()
{
    clearTreeCanonInterner();
    const vector<edgeSpec> triangle = cycleEdges(3);
    molGraph first = makeGraph(vector<string>(3, "C"), triangle);
    molGraph second = makeGraph(
        vector<string>(3, "C"),
        triangle,
        shuffledPermutation(3, 101),
        true
    );
    molGraph labelled = makeGraph(vector<string>{"C", "C", "N"}, triangle);

    const graphHash firstHash(first, true);
    const graphHash secondHash(second, true);
    const graphHash labelledHash(labelled, true);
    const size_t storedHash = std::hash<graphHash>{}(firstHash);
    require(
        firstHash.cyclicHash.canonicalCode.empty(),
        "cyclic graph hash did not begin with a lazy exact code"
    );
    firstHash.prepareForSharing();
    require(
        !firstHash.cyclicHash.canonicalCode.empty(),
        "shared cyclic graph hash did not materialise its exact code"
    );

    first = molGraph{};
    second = molGraph{};
    labelled = molGraph{};
    targetMolecule = molGraph{};
    univEdgeList.clear();
    clearTreeCanonInterner();

    require(firstHash == secondHash, "cached equal graph hashes changed with source");
    require(!(firstHash == labelledHash), "cached unequal graph hashes lost labels");
    require(
        std::hash<graphHash>{}(firstHash) == storedHash,
        "cached graph hash changed after source/interner destruction"
    );

    std::unordered_set<graphHash> cachedKeys;
    require(cachedKeys.insert(firstHash).second, "first cached key was rejected");
    require(
        !cachedKeys.insert(secondHash).second,
        "isomorphic cached key was inserted twice"
    );
    require(
        cachedKeys.insert(labelledHash).second,
        "non-isomorphic cached key was merged"
    );
}

void configureCanoniseMaskDomain(size_t edgeCount)
{
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(edgeCount);
    std::construct_at(std::addressof(allEdges));
}

EdgeMask componentMask(int firstVertex, int vertexCount)
{
    EdgeMask result;
    const int lastVertex = firstVertex + vertexCount;
    const vector<edgeL> &edgeList = searchUniverseEdgeList();
    for (size_t edge = 0; edge < edgeList.size(); edge++)
    {
        if (
            edgeList[edge].a >= firstVertex &&
            edgeList[edge].a < lastVertex &&
            edgeList[edge].b >= firstVertex &&
            edgeList[edge].b < lastVertex
        ) result.set(edge);
    }
    return result;
}

void testSharedCanonicalRegistryAdversarialOrdering()
{
    sharedCanonicalIdRegistry registry(37);
    constexpr int workerCount = 8;
    std::barrier start(workerCount);
    std::array<int, workerCount> bipartiteIds{};
    std::array<int, workerCount> prismIds{};
    std::atomic<int> invariantFailures{0};
    std::vector<std::thread> workers;
    workers.reserve(workerCount);

    for (int worker = 0; worker < workerCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            bitsetHashTable.clear();
            graphHashMap.clear();
            clearTreeCanonInterner();
            targetMolecule = makeSharedRegistryCollisionGraph();
            univEdgeList = targetMolecule.writeEdgeList();
            configureCanoniseMaskDomain(univEdgeList.size());
            prepareCanonicalisationGraph(targetMolecule, univEdgeList);
            sharedCanonicalRegistry = &registry;

            {
                // Alternate exact occurrences so equivalent candidates also
                // arrive with different serialized representative masks.
                EdgeMask bipartite = componentMask(
                    worker % 2 == 0 ? 0 : 6,
                    6
                );
                EdgeMask prism = componentMask(
                    worker % 2 == 0 ? 12 : 18,
                    6
                );
                bool bipartiteIsCyclic = false;
                const flatCanonGraph &bipartiteGraph =
                    canonicalisationGraphScratch.build(
                        targetMolecule,
                        univEdgeList,
                        bipartite,
                        bipartiteIsCyclic
                    );
                const size_t bipartiteInvariant =
                    sharedCanonicalInvariantHash(bipartiteGraph);
                bool prismIsCyclic = false;
                const flatCanonGraph &prismGraph =
                    canonicalisationGraphScratch.build(
                        targetMolecule,
                        univEdgeList,
                        prism,
                        prismIsCyclic
                    );
                const size_t prismInvariant =
                    sharedCanonicalInvariantHash(prismGraph);
                if (
                    !bipartiteIsCyclic || !prismIsCyclic ||
                    bipartiteInvariant != prismInvariant
                )
                {
                    invariantFailures.fetch_add(
                        1,
                        std::memory_order_relaxed
                    );
                }

                start.arrive_and_wait();
                if (worker % 2 == 0)
                {
                    bipartiteIds[worker] = canonise(bipartite);
                    prismIds[worker] = canonise(prism);
                }
                else
                {
                    prismIds[worker] = canonise(prism);
                    bipartiteIds[worker] = canonise(bipartite);
                }
            }

            sharedCanonicalRegistry = nullptr;
            bitsetHashTable.clear();
            graphHashMap.clear();
            clearTreeCanonInterner();
        });
    }
    for (std::thread &worker : workers) worker.join();

    require(
        invariantFailures.load(std::memory_order_relaxed) == 0,
        "registry collision fixtures did not share one invariant bucket"
    );
    require(
        std::all_of(
            bipartiteIds.begin(),
            bipartiteIds.end(),
            [&](int id) {return id == bipartiteIds.front();}
        ),
        "equivalent K3,3 candidates received different shared IDs"
    );
    require(
        std::all_of(
            prismIds.begin(),
            prismIds.end(),
            [&](int id) {return id == prismIds.front();}
        ),
        "equivalent prism candidates received different shared IDs"
    );
    require(
        bipartiteIds.front() != prismIds.front(),
        "non-isomorphic graphs sharing an invariant bucket were merged"
    );
    require(
        bipartiteIds.front() >= 37 && prismIds.front() >= 37,
        "shared registry allocated below its first unused ID"
    );
}

void testSharedCanonicalRegistryWithDivergentTreeInterners()
{
    sharedCanonicalRegistry = nullptr;
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    targetMolecule = makeSharedRegistryTreeInternerGraph();
    univEdgeList = targetMolecule.writeEdgeList();
    configureCanoniseMaskDomain(univEdgeList.size());
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);

    EdgeMask producerSeed = componentMask(36, 2);
    const int producerSeedId = canonise(producerSeed);
    require(
        producerSeedId == 0,
        "producer labelled-edge seed did not receive the first ID"
    );
    EdgeMask cyclicProducerSeed = componentMask(38, 3);
    const int cyclicProducerSeedId = canonise(cyclicProducerSeed);
    require(
        cyclicProducerSeedId == 1,
        "producer cyclic seed did not receive the second ID"
    );
    require(
        graphHashMap.size() == 2,
        "producer did not create exactly two canonical seed classes"
    );

    freezeGraphHashSeed(graphHashMap);
    const auto seedGraphHashes = graphHashMap;
    const auto seedAtomInterner = treeCanonAtomInterner;
    const auto seedLeafInterner = treeCanonLeafInterner;
    const auto seedTreeInterner = treeCanonInterner;
    const molGraph sharedMolecule = std::move(targetMolecule);
    const vector<edgeL> sharedEdgeList = std::move(univEdgeList);
    vector<uint64_t> frozenCyclicSeedCode;
    for (const auto &entry : seedGraphHashes)
    {
        if (entry.second.first == cyclicProducerSeedId)
            frozenCyclicSeedCode = entry.first.cyclicHash.canonicalCode;
    }
    require(
        !frozenCyclicSeedCode.empty(),
        "cyclic graph-hash seed was not frozen before publication"
    );
    sharedCanonicalIdRegistry registry(graphHashMap.size());

    constexpr int workerCount = 8;
    std::barrier start(workerCount);
    std::array<int, workerCount> forkedTreeIds{};
    std::array<int, workerCount> bentArmTreeIds{};
    std::array<int, workerCount> pendantCycleIds{};
    std::array<int, workerCount> producerSeedIds{};
    std::array<int, workerCount> cyclicProducerSeedIds{};
    std::array<treeCanonNodeId, workerCount> singleBondSignatureIds{};
    std::array<size_t, workerCount> graphDeltaSizes{};
    std::array<size_t, workerCount> treeDeltaSizes{};
    std::array<size_t, workerCount> leafDeltaSizes{};
    std::atomic<int> invariantFailures{0};
    std::vector<std::thread> workers;
    workers.reserve(workerCount);

    for (int worker = 0; worker < workerCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            bitsetHashTable.clear();
            clearGraphHashDelta();
            clearTreeCanonInterner();
            sharedTargetMolecule = &sharedMolecule;
            sharedUniverseEdgeList = &sharedEdgeList;
            configureCanoniseMaskDomain(sharedEdgeList.size());
            bindTreeCanonInternerSeed(
                seedAtomInterner,
                seedLeafInterner,
                seedTreeInterner
            );
            bindGraphHashSeed(seedGraphHashes);
            prepareCanonicalisationGraph(
                searchTargetMolecule(),
                searchUniverseEdgeList()
            );
            if (
                std::addressof(searchTargetMolecule()) !=
                    std::addressof(sharedMolecule) ||
                searchUniverseEdgeList().data() != sharedEdgeList.data() ||
                sharedTreeCanonAtomInterner != &seedAtomInterner ||
                sharedTreeCanonLeafInterner != &seedLeafInterner ||
                sharedTreeCanonInterner != &seedTreeInterner ||
                sharedGraphHashSeed != &seedGraphHashes ||
                !graphHashMap.empty() || !treeCanonInterner.empty() ||
                !treeCanonLeafInternerDelta.empty()
            )
            {
                invariantFailures.fetch_add(1, std::memory_order_relaxed);
            }

            // Allocate two exact rooted-tree signatures in opposite orders.
            // The single-bond form is used by both tree fixtures and by the
            // pendant chain, so its numeric TLS ID must genuinely diverge.
            const treeCanonAtomId carbon = internTreeCanonAtom("C");
            const treeCanonNodeId leaf = internTreeCanonNode(
                carbon,
                vector<treeCanonChild>{}
            );
            const auto internBondSignature = [&](uint16_t bondType)
            {
                return internTreeCanonNode(
                    carbon,
                    vector<treeCanonChild>{{bondType, leaf}}
                );
            };
            if (worker % 2 == 0)
            {
                singleBondSignatureIds[worker] = internBondSignature(1);
                static_cast<void>(internBondSignature(2));
            }
            else
            {
                static_cast<void>(internBondSignature(2));
                singleBondSignatureIds[worker] = internBondSignature(1);
            }

            sharedCanonicalRegistry = &registry;
            {
                EdgeMask forkedTree = componentMask(
                    worker % 2 == 0 ? 0 : 6,
                    6
                );
                EdgeMask bentArmTree = componentMask(
                    worker % 2 == 0 ? 12 : 18,
                    6
                );
                EdgeMask pendantCycle = componentMask(
                    worker % 2 == 0 ? 24 : 30,
                    6
                );
                EdgeMask seededClass = componentMask(36, 2);
                EdgeMask cyclicSeededClass = componentMask(38, 3);

                bool forkedIsCyclic = false;
                const flatCanonGraph &forkedGraph =
                    canonicalisationGraphScratch.build(
                        searchTargetMolecule(),
                        searchUniverseEdgeList(),
                        forkedTree,
                        forkedIsCyclic
                    );
                const size_t forkedInvariant =
                    sharedCanonicalInvariantHash(forkedGraph);
                bool bentArmIsCyclic = false;
                const flatCanonGraph &bentArmGraph =
                    canonicalisationGraphScratch.build(
                        searchTargetMolecule(),
                        searchUniverseEdgeList(),
                        bentArmTree,
                        bentArmIsCyclic
                    );
                const size_t bentArmInvariant =
                    sharedCanonicalInvariantHash(bentArmGraph);
                if (
                    forkedIsCyclic || bentArmIsCyclic ||
                    forkedInvariant != bentArmInvariant
                )
                {
                    invariantFailures.fetch_add(
                        1,
                        std::memory_order_relaxed
                    );
                }

                start.arrive_and_wait();
                if (worker % 2 == 0)
                {
                    forkedTreeIds[worker] = canonise(forkedTree);
                    bentArmTreeIds[worker] = canonise(bentArmTree);
                }
                else
                {
                    bentArmTreeIds[worker] = canonise(bentArmTree);
                    forkedTreeIds[worker] = canonise(forkedTree);
                }
                pendantCycleIds[worker] = canonise(pendantCycle);
                producerSeedIds[worker] = canonise(seededClass);
                cyclicProducerSeedIds[worker] = canonise(cyclicSeededClass);
            }

            graphDeltaSizes[worker] = graphHashMap.size();
            treeDeltaSizes[worker] = treeCanonInterner.size();
            leafDeltaSizes[worker] = treeCanonLeafInternerDelta.size();

            sharedCanonicalRegistry = nullptr;
            bitsetHashTable.clear();
            clearGraphHashDelta();
            clearTreeCanonInterner();
            sharedTargetMolecule = nullptr;
            sharedUniverseEdgeList = nullptr;
        });
    }
    for (std::thread &worker : workers) worker.join();

    require(
        singleBondSignatureIds[0] != singleBondSignatureIds[1],
        "test workers did not diverge their TLS tree-signature IDs"
    );
    for (int worker = 0; worker < workerCount; ++worker)
    {
        const treeCanonNodeId expected =
            singleBondSignatureIds[worker % 2];
        require(
            singleBondSignatureIds[worker] == expected,
            "same-order workers allocated inconsistent TLS signatures"
        );
    }
    require(
        invariantFailures.load(std::memory_order_relaxed) == 0,
        "shared seed binding or invariant bucket was inconsistent"
    );
    require(
        std::all_of(
            graphDeltaSizes.begin(),
            graphDeltaSizes.end(),
            [](size_t size) {return size != 0;}
        ),
        "workers did not retain graph hashes in local deltas"
    );
    require(
        std::all_of(
            treeDeltaSizes.begin(),
            treeDeltaSizes.end(),
            [](size_t size) {return size != 0;}
        ),
        "workers did not retain tree signatures in local deltas"
    );
    require(
        std::all_of(
            leafDeltaSizes.begin(),
            leafDeltaSizes.end(),
            [](size_t size) {return size != 0;}
        ),
        "workers did not retain leaves in local deltas"
    );
    require(
        std::all_of(
            forkedTreeIds.begin(),
            forkedTreeIds.end(),
            [&](int id) {return id == forkedTreeIds.front();}
        ),
        "equivalent dynamic forked trees received different shared IDs"
    );
    require(
        std::all_of(
            bentArmTreeIds.begin(),
            bentArmTreeIds.end(),
            [&](int id) {return id == bentArmTreeIds.front();}
        ),
        "equivalent dynamic tree candidates received different shared IDs"
    );
    require(
        bentArmTreeIds.front() != forkedTreeIds.front() &&
            bentArmTreeIds.front() != producerSeedId,
        "non-isomorphic dynamic tree class aliased a producer seed"
    );
    require(
        forkedTreeIds.front() >= 2 && bentArmTreeIds.front() >= 2 &&
            pendantCycleIds.front() >= 2,
        "dynamic class reused the producer's canonical ID range"
    );
    require(
        std::all_of(
            pendantCycleIds.begin(),
            pendantCycleIds.end(),
            [&](int id) {return id == pendantCycleIds.front();}
        ),
        "dynamic pendant-cycle class changed across divergent TLS interners"
    );
    require(
        pendantCycleIds.front() != forkedTreeIds.front() &&
            pendantCycleIds.front() != bentArmTreeIds.front() &&
            pendantCycleIds.front() != producerSeedId,
        "dynamic pendant-cycle class aliased another canonical class"
    );
    require(
        std::all_of(
            producerSeedIds.begin(),
            producerSeedIds.end(),
            [&](int id) {return id == producerSeedId;}
        ),
        "producer-seeded class changed after dynamic ID allocation"
    );
    require(
        std::all_of(
            cyclicProducerSeedIds.begin(),
            cyclicProducerSeedIds.end(),
            [&](int id) {return id == cyclicProducerSeedId;}
        ),
        "frozen cyclic seed changed during concurrent lookup"
    );
    bool cyclicSeedUnchanged = false;
    bool treeSeedUnchanged = false;
    for (const auto &entry : seedGraphHashes)
    {
        if (entry.second == pii{producerSeedId, 1})
            treeSeedUnchanged = true;
        if (
            entry.second == pii{cyclicProducerSeedId, 1} &&
            entry.first.cyclicHash.canonicalCode == frozenCyclicSeedCode
        ) cyclicSeedUnchanged = true;
    }
    require(
        seedGraphHashes.size() == 2 && treeSeedUnchanged &&
            cyclicSeedUnchanged && seedAtomInterner.size() == 4 &&
            !seedLeafInterner.empty() &&
            !seedTreeInterner.empty(),
        "workers mutated the shared producer seed"
    );

    sharedCanonicalRegistry = nullptr;
    bitsetHashTable.clear();
    clearGraphHashDelta();
    clearTreeCanonInterner();
    sharedTargetMolecule = nullptr;
    sharedUniverseEdgeList = nullptr;
}

void testCanonicalSeedOverlayWithoutRegistry()
{
    sharedCanonicalRegistry = nullptr;
    bitsetHashTable.clear();
    clearGraphHashDelta();
    clearTreeCanonInterner();
    targetMolecule = makeSharedRegistryTreeInternerGraph();
    univEdgeList = targetMolecule.writeEdgeList();
    configureCanoniseMaskDomain(univEdgeList.size());
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);

    int producerSeedId;
    {
        EdgeMask producerSeed = componentMask(36, 2);
        producerSeedId = canonise(producerSeed);
    }
    require(producerSeedId == 0, "overlay seed did not receive ID zero");
    freezeGraphHashSeed(graphHashMap);

    const auto seedGraphHashes = graphHashMap;
    const auto seedAtomInterner = treeCanonAtomInterner;
    const auto seedLeafInterner = treeCanonLeafInterner;
    const auto seedTreeInterner = treeCanonInterner;
    const molGraph sharedMolecule = std::move(targetMolecule);
    const vector<edgeL> sharedEdgeList = std::move(univEdgeList);

    bitsetHashTable.clear();
    clearGraphHashDelta();
    clearTreeCanonInterner();
    sharedTargetMolecule = &sharedMolecule;
    sharedUniverseEdgeList = &sharedEdgeList;
    bindTreeCanonInternerSeed(
        seedAtomInterner,
        seedLeafInterner,
        seedTreeInterner
    );
    bindGraphHashSeed(seedGraphHashes);
    prepareCanonicalisationGraph(
        searchTargetMolecule(),
        searchUniverseEdgeList()
    );

    {
        EdgeMask producerSeed = componentMask(36, 2);
        EdgeMask firstFork = componentMask(0, 6);
        EdgeMask secondFork = componentMask(6, 6);
        EdgeMask bentArm = componentMask(12, 6);
        require(
            canonise(producerSeed) == producerSeedId,
            "overlay lookup changed a producer canonical ID"
        );
        const int forkId = canonise(firstFork);
        require(
            forkId == 1 && canonise(secondFork) == forkId,
            "overlay did not allocate/reuse the first local canonical ID"
        );
        require(
            canonise(bentArm) == 2,
            "overlay local canonical IDs overlapped the seed range"
        );
    }
    require(
        graphHashMap.size() == 2 && graphHashClassCount() == 3 &&
            seedGraphHashes.size() == 1 &&
            seedGraphHashes.begin()->second == pii{producerSeedId, 1},
        "no-registry overlay mutated or miscounted its shared seed"
    );

    bitsetHashTable.clear();
    clearGraphHashDelta();
    clearTreeCanonInterner();
    sharedTargetMolecule = nullptr;
    sharedUniverseEdgeList = nullptr;
}

void testFlatCanoniseMaskPath()
{
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();

    const vector<string> labels{
        "C", "C", "N",
        "N", "C", "C",
        "C", "C", "O"
    };
    const vector<edgeSpec> edges{
        {0, 1, 1}, {1, 2, 2}, {2, 0, 1},
        {3, 4, 2}, {4, 5, 1}, {5, 3, 1},
        {6, 7, 1}, {7, 8, 2}, {8, 6, 1}
    };
    targetMolecule = makeGraph(labels, edges);
    univEdgeList = targetMolecule.writeEdgeList();
    configureCanoniseMaskDomain(univEdgeList.size());
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);

    {
        EdgeMask first = componentMask(0, 3);
        EdgeMask isomorphic = componentMask(3, 3);
        EdgeMask different = componentMask(6, 3);
        require(first.count() == 3, "first flat canonical mask is incomplete");
        require(
            canonise(first) == canonise(isomorphic),
            "flat cyclic canonicalisation is vertex-order dependent"
        );
        require(
            canonise(first) != canonise(different),
            "flat cyclic canonicalisation lost an atom label"
        );
    }

    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();

    targetMolecule = makeGraph(
        {"C", "C", "N", "O", "N", "C", "C", "O"},
        {
            {0, 1, 1}, {1, 2, 2}, {2, 0, 1}, {2, 3, 1},
            {5, 6, 1}, {6, 4, 2}, {4, 5, 1}, {4, 7, 1}
        }
    );
    univEdgeList = targetMolecule.writeEdgeList();
    configureCanoniseMaskDomain(univEdgeList.size());
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
    {
        EdgeMask first = componentMask(0, 4);
        EdgeMask intervening = componentMask(0, 3);
        EdgeMask isomorphic = componentMask(4, 4);
        const int firstId = canonise(first);
        (void)canonise(intervening);
        require(
            firstId == canonise(isomorphic),
            "flat cyclic peeling lost a pendant-tree canonical form"
        );
    }

    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();

    targetMolecule = makeGraph(
        {
            "C", "C", "N", "O", "S",
            "N", "C", "C", "S", "O"
        },
        {
            {0, 1, 1}, {1, 2, 2}, {2, 0, 1}, {3, 4, 1},
            {6, 7, 1}, {7, 5, 2}, {5, 6, 1}, {8, 9, 1}
        }
    );
    univEdgeList = targetMolecule.writeEdgeList();
    configureCanoniseMaskDomain(univEdgeList.size());
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
    {
        EdgeMask first = componentMask(0, 5);
        EdgeMask intervening = componentMask(3, 2);
        EdgeMask isomorphic = componentMask(5, 5);
        const int firstId = canonise(first);
        (void)canonise(intervening);
        require(
            firstId == canonise(isomorphic),
            "flat disconnected fallback retained reusable CSR storage"
        );
    }

    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();

    targetMolecule = makeGraph(
        {"C", "N"},
        {{0, 1, 1}, {0, 1, 2}}
    );
    univEdgeList = targetMolecule.writeEdgeList();
    configureCanoniseMaskDomain(univEdgeList.size());
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
    int parallelBondId;
    {
        EdgeMask parallelBonds;
        parallelBonds.set();
        parallelBondId = canonise(parallelBonds);
    }
    bitsetHashTable.clear();
    targetMolecule = makeGraph(
        {"N", "C"},
        {{0, 1, 1}, {0, 1, 3}}
    );
    univEdgeList = targetMolecule.writeEdgeList();
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
    {
        EdgeMask changedParallelBond;
        changedParallelBond.set();
        require(
            parallelBondId != canonise(changedParallelBond),
            "flat cyclic canonicalisation lost a parallel bond label"
        );
    }

    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    targetMolecule = makeGraph(
        {"C", "X", "C"},
        {{0, 1, 1}, {1, 2, 1}}
    );
    univEdgeList = targetMolecule.writeEdgeList();
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
    {
        EdgeMask legacyPath;
        legacyPath.set();
        bool isCyclic = false;
        const flatCanonGraph &flat = canonicalisationGraphScratch.build(
            targetMolecule,
            univEdgeList,
            legacyPath,
            isCyclic
        );
        const graphHash hash(flat, isCyclic);
        require(flat.hasLegacyX, "flat canonical graph lost the X sentinel");
        require(
            hash.treeHash.empty() && !hash.cyclicHash.empty(),
            "flat X-sentinel graph did not use whole-graph fallback"
        );
    }

    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();

    vector<string> pathLabels(67, "C");
    vector<edgeSpec> pathEdges;
    pathEdges.reserve(66);
    for (int edge = 0; edge < 66; edge++)
        pathEdges.emplace_back(edge, edge + 1, 1);
    targetMolecule = makeGraph(pathLabels, pathEdges);
    univEdgeList = targetMolecule.writeEdgeList();
    configureCanoniseMaskDomain(univEdgeList.size());
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);

    {
        EdgeMask lowEdge;
        EdgeMask highEdge;
        lowEdge.set(0);
        highEdge.set(65);
        require(
            canonise(lowEdge) == canonise(highEdge),
            "flat canonicalisation missed a selected wide-mask bit"
        );
    }

    bitsetHashTable.clear();
    graphHashMap.clear();
    configureCanoniseMaskDomain(64);
}

int main()
{
    testAgainstExactMatcher();
    testLargeAttachedTrees();
    testSharedCanonicalRegistryAdversarialOrdering();
    testSharedCanonicalRegistryWithDivergentTreeInterners();
    testCanonicalSeedOverlayWithoutRegistry();
    testFlatCanoniseMaskPath();
    testCachedGraphHashIsSelfContained();
    graphHashMap.clear();
    clearTreeCanonInterner();
    return 0;
}
