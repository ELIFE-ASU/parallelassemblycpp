#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using treeCanonNodeId = std::uint64_t;
using treeCanonAtomId = std::uint64_t;

/** One directed edge in the reusable canonicalisation CSR input. */
struct flatCanonAdjacentEdge
{
    std::uint32_t neighbour = 0;
    std::uint16_t bondType = 0;
};

/**
 * @brief Flat edge-induced graph passed directly to the canonicalisers.
 *
 * The vectors are reusable miss-path scratch. Canonical forms must therefore
 * copy any representation retained after the current canonicalisation call.
 */
struct flatCanonGraph
{
    std::vector<treeCanonAtomId> labels;
    std::vector<std::size_t> adjacencyOffsets;
    std::vector<flatCanonAdjacentEdge> adjacency;
    std::size_t edgeCount = 0;
    bool hasLegacyX = false;
    bool hasPendantVertex = false;

    [[nodiscard]] std::span<const flatCanonAdjacentEdge> neighbours(
        std::size_t vertex
    ) const noexcept
    {
        return std::span<const flatCanonAdjacentEdge>(adjacency).subspan(
            adjacencyOffsets[vertex],
            adjacencyOffsets[vertex + 1] - adjacencyOffsets[vertex]
        );
    }
};

/**
 * @brief Canonical identity of an unrooted, labelled tree.
 *
 * A tree has either one centroid, represented by first alone, or two
 * centroids, represented by the sorted first/second pair and their connecting
 * bond. Zero is reserved for graphs without a tree canonical form.
 */
struct treeCanonForm
{
    treeCanonNodeId first = 0;
    treeCanonNodeId second = 0;
    std::uint16_t centralBond = 0;

    bool empty() const {return first == 0;}
    bool operator==(const treeCanonForm &) const = default;
};

/**
 * @brief One labelled edge to an already-canonical child subtree.
 */
struct treeCanonChild
{
    std::uint16_t bondType;
    treeCanonNodeId subtree;

    bool operator==(const treeCanonChild &) const = default;
};

struct treeCanonWorkspace
{
    std::vector<int> remainingDegree;
    std::vector<unsigned char> removed;
    std::vector<std::vector<treeCanonChild>> children;
    std::vector<int> leaves;
    std::vector<int> nextLeaves;

    void begin(std::size_t nodeCount)
    {
        if (remainingDegree.size() < nodeCount)
            remainingDegree.resize(nodeCount);
        if (removed.size() < nodeCount) removed.resize(nodeCount);
        std::fill_n(removed.begin(), nodeCount, 0);
        if (children.size() < nodeCount) children.resize(nodeCount);
        for (std::size_t node = 0; node < nodeCount; node++)
            children[node].clear();
        leaves.clear();
        nextLeaves.clear();
        if (leaves.capacity() < nodeCount) leaves.reserve(nodeCount);
        if (nextLeaves.capacity() < nodeCount) nextLeaves.reserve(nodeCount);
    }
};

// Canonical maps and the surrounding search are already single-threaded.
inline ASSEMBLYCPP_SEARCH_LOCAL treeCanonWorkspace treeCanonScratch;

/**
 * @brief Exact structural key for a rooted subtree.
 */
struct treeCanonSignature
{
    treeCanonAtomId atomType;
    std::vector<treeCanonChild> children;

    bool operator==(const treeCanonSignature &) const = default;
};

struct treeCanonSignatureHash
{
    static void combine(std::size_t &seed, std::size_t value)
    {
        seed ^= value + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL)
            + (seed << 6) + (seed >> 2);
    }

    std::size_t operator()(const treeCanonSignature &signature) const
    {
        std::size_t result = std::hash<treeCanonAtomId>{}(signature.atomType);
        for (const treeCanonChild &child : signature.children)
        {
            combine(result, std::hash<std::uint16_t>{}(child.bondType));
            combine(result, std::hash<treeCanonNodeId>{}(child.subtree));
        }
        return result;
    }
};

/**
 * IDs are opaque and stable within one interner generation. Production
 * discards graphHashMap before resetting these interners between calculations;
 * forms from different generations must not be mixed.
 */
inline ASSEMBLYCPP_SEARCH_LOCAL std::unordered_map<
    std::string,
    treeCanonAtomId
> treeCanonAtomInterner;
inline ASSEMBLYCPP_SEARCH_LOCAL std::vector<treeCanonNodeId> treeCanonLeafInterner;
inline std::unordered_map<
    treeCanonSignature,
    treeCanonNodeId,
    treeCanonSignatureHash
> ASSEMBLYCPP_SEARCH_LOCAL treeCanonInterner;
inline ASSEMBLYCPP_SEARCH_LOCAL std::uint64_t treeCanonInternerGeneration = 1;

void clearTreeCanonInterner()
{
    treeCanonInterner.clear();
    treeCanonAtomInterner.clear();
    treeCanonLeafInterner.clear();
    if (++treeCanonInternerGeneration == 0) treeCanonInternerGeneration = 1;
}

treeCanonAtomId internTreeCanonAtom(const std::string &atomType)
{
    const treeCanonAtomId nextId =
        static_cast<treeCanonAtomId>(treeCanonAtomInterner.size()) + 1;
    return treeCanonAtomInterner.try_emplace(atomType, nextId).first->second;
}

treeCanonNodeId internTreeCanonSignature(treeCanonSignature signature)
{
    const treeCanonNodeId nextId =
        static_cast<treeCanonNodeId>(treeCanonInterner.size()) + 1;
    return treeCanonInterner.try_emplace(
        std::move(signature),
        nextId
    ).first->second;
}

treeCanonNodeId internTreeCanonNode(
    treeCanonAtomId atomType,
    std::vector<treeCanonChild> children
)
{
    std::sort(
        children.begin(),
        children.end(),
        [](const treeCanonChild &left, const treeCanonChild &right)
        {
            if (left.bondType != right.bondType)
                return left.bondType < right.bondType;
            return left.subtree < right.subtree;
        }
    );
    if (children.empty())
    {
        if (treeCanonLeafInterner.size() <= atomType)
            treeCanonLeafInterner.resize(atomType + 1, 0);
        treeCanonNodeId &leaf = treeCanonLeafInterner[atomType];
        if (leaf == 0)
            leaf = internTreeCanonSignature({atomType, {}});
        return leaf;
    }
    return internTreeCanonSignature({atomType, std::move(children)});
}

inline treeCanonNodeId internTreeCanonNode(
    molGraph &mg,
    int node,
    std::vector<treeCanonChild> children
)
{
    return internTreeCanonNode(
        internTreeCanonAtom(mg.mg[node].type),
        std::move(children)
    );
}

[[nodiscard]] inline std::size_t canonGraphVertexCount(
    const molGraph &graph
) noexcept
{
    return graph.mg.size();
}

[[nodiscard]] inline std::size_t canonGraphVertexCount(
    const flatCanonGraph &graph
) noexcept
{
    return graph.labels.size();
}

[[nodiscard]] inline bool canonGraphHasEdgeCount(
    const molGraph &graph,
    std::size_t expected
) noexcept
{
    return graph.totalBonds >= 0 &&
        static_cast<std::size_t>(graph.totalBonds) == expected;
}

[[nodiscard]] inline bool canonGraphHasEdgeCount(
    const flatCanonGraph &graph,
    std::size_t expected
) noexcept
{
    return graph.edgeCount == expected;
}

[[nodiscard]] inline bool canonGraphHasLegacyX(const molGraph &graph)
{
    return std::any_of(
        graph.mg.begin(),
        graph.mg.end(),
        [](const atom &vertex) {return vertex.type == "X";}
    );
}

[[nodiscard]] inline bool canonGraphHasLegacyX(
    const flatCanonGraph &graph
) noexcept
{
    return graph.hasLegacyX;
}

[[nodiscard]] inline bool canonGraphHasPendantVertex(const molGraph &graph)
{
    return std::any_of(
        graph.mg.begin(),
        graph.mg.end(),
        [](const atom &vertex) {return vertex.list.size() <= 1;}
    );
}

[[nodiscard]] inline bool canonGraphHasPendantVertex(
    const flatCanonGraph &graph
) noexcept
{
    return graph.hasPendantVertex;
}

[[nodiscard]] inline treeCanonAtomId canonGraphAtomType(
    const molGraph &graph,
    std::size_t vertex
)
{
    return internTreeCanonAtom(graph.mg[vertex].type);
}

[[nodiscard]] inline treeCanonAtomId canonGraphAtomType(
    const flatCanonGraph &graph,
    std::size_t vertex
) noexcept
{
    return graph.labels[vertex];
}

[[nodiscard]] inline const std::vector<bond> &canonGraphNeighbours(
    const molGraph &graph,
    std::size_t vertex
) noexcept
{
    return graph.mg[vertex].list;
}

[[nodiscard]] inline std::span<const flatCanonAdjacentEdge> canonGraphNeighbours(
    const flatCanonGraph &graph,
    std::size_t vertex
) noexcept
{
    return graph.neighbours(vertex);
}

[[nodiscard]] inline int canonGraphNeighbour(const bond &edge) noexcept
{
    return edge.n;
}

[[nodiscard]] inline int canonGraphNeighbour(
    const flatCanonAdjacentEdge &edge
) noexcept
{
    return static_cast<int>(edge.neighbour);
}

[[nodiscard]] inline std::uint16_t canonGraphBondType(const bond &edge) noexcept
{
    return static_cast<std::uint16_t>(edge.type);
}

[[nodiscard]] inline std::uint16_t canonGraphBondType(
    const flatCanonAdjacentEdge &edge
) noexcept
{
    return edge.bondType;
}

/**
 * @brief Iteratively canonicalise an acyclic molecular graph.
 *
 * AHU leaf peeling processes a complete layer at a time. Each removed node is
 * represented by an exact interned ID built from its atom label and sorted
 * (bond-label, child-ID) multiset. The one or two unremoved nodes are the tree
 * centroids, so no recursive traversal or subtree-string concatenation is
 * required.
 *
 * @param mg Acyclic, connected molGraph
 * @param n Retained for API compatibility; centroid selection is root-free.
 */
template<typename Graph>
treeCanonForm centroidTreeCanonImpl(const Graph &graph, int n)
{
    (void)n;
    const std::size_t nodeCount = canonGraphVertexCount(graph);
    if (nodeCount == 0) return {};
    if (!canonGraphHasEdgeCount(graph, nodeCount - 1)) return {};
    if (canonGraphHasLegacyX(graph)) return {};

    treeCanonScratch.begin(nodeCount);
    auto &remainingDegree = treeCanonScratch.remainingDegree;
    auto &removed = treeCanonScratch.removed;
    auto &children = treeCanonScratch.children;
    auto &leaves = treeCanonScratch.leaves;
    auto &nextLeaves = treeCanonScratch.nextLeaves;
    std::size_t remaining = nodeCount;

    for (std::size_t node = 0; node < nodeCount; node++)
    {
        remainingDegree[node] = static_cast<int>(
            canonGraphNeighbours(graph, node).size()
        );
        if (remainingDegree[node] <= 1)
            leaves.push_back(static_cast<int>(node));
    }

    while (remaining > 2)
    {
        // A connected acyclic graph always has a leaf. Returning no tree form
        // safely routes malformed input through exact whole-graph canonicalisation.
        if (leaves.empty()) return {};
        remaining -= leaves.size();
        for (const int leaf : leaves) removed[leaf] = 1;

        nextLeaves.clear();
        for (const int leaf : leaves)
        {
            const treeCanonNodeId subtree = internTreeCanonNode(
                canonGraphAtomType(graph, static_cast<std::size_t>(leaf)),
                std::move(children[leaf])
            );
            for (const auto &edge : canonGraphNeighbours(graph, leaf))
            {
                const int neighbour = canonGraphNeighbour(edge);
                if (
                    neighbour < 0 ||
                    static_cast<std::size_t>(neighbour) >= nodeCount ||
                    removed[neighbour]
                ) continue;

                children[neighbour].push_back({
                    canonGraphBondType(edge),
                    subtree
                });
                remainingDegree[neighbour]--;
                if (remainingDegree[neighbour] == 1)
                    nextLeaves.push_back(neighbour);
                break;
            }
        }
        leaves.swap(nextLeaves);
    }

    std::array<int, 2> centroids{};
    std::size_t centroidCount = 0;
    for (std::size_t node = 0; node < nodeCount; node++)
    {
        if (removed[node]) continue;
        if (centroidCount == centroids.size()) return {};
        centroids[centroidCount++] = static_cast<int>(node);
    }
    if (centroidCount == 0) return {};

    const treeCanonNodeId first = internTreeCanonNode(
        canonGraphAtomType(graph, static_cast<std::size_t>(centroids[0])),
        std::move(children[centroids[0]])
    );
    if (centroidCount == 1) return {first, 0, 0};

    std::uint16_t centralBond = 0;
    bool centralBondFound = false;
    for (const auto &edge : canonGraphNeighbours(graph, centroids[0]))
    {
        if (canonGraphNeighbour(edge) == centroids[1])
        {
            centralBond = canonGraphBondType(edge);
            centralBondFound = true;
            break;
        }
    }
    if (!centralBondFound) return {};

    const treeCanonNodeId second = internTreeCanonNode(
        canonGraphAtomType(graph, static_cast<std::size_t>(centroids[1])),
        std::move(children[centroids[1]])
    );
    if (first <= second) return {first, second, centralBond};
    return {second, first, centralBond};
}

inline treeCanonForm centroidTreeCanon(molGraph &mg, int n)
{
    return centroidTreeCanonImpl(mg, n);
}

inline treeCanonForm centroidTreeCanon(const flatCanonGraph &graph, int n)
{
    return centroidTreeCanonImpl(graph, n);
}
