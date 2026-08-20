#include <algorithm>
#include <array>
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

int main()
{
    testAgainstExactMatcher();
    testLargeAttachedTrees();
    testCachedGraphHashIsSelfContained();
    graphHashMap.clear();
    clearTreeCanonInterner();
    return 0;
}
