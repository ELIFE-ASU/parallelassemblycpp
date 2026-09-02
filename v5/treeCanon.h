#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
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

// Each OpenMP worker owns this scratch through ASSEMBLYCPP_SEARCH_LOCAL.
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
// Parallel workers layer small mutable interners over the producer's frozen
// generation. The seed pointers are borrowed from SearchContext for exactly
// the lifetime of the parallel region.
inline ASSEMBLYCPP_SEARCH_LOCAL const decltype(treeCanonAtomInterner)
    *sharedTreeCanonAtomInterner = nullptr;
inline ASSEMBLYCPP_SEARCH_LOCAL const decltype(treeCanonLeafInterner)
    *sharedTreeCanonLeafInterner = nullptr;
inline ASSEMBLYCPP_SEARCH_LOCAL const decltype(treeCanonInterner)
    *sharedTreeCanonInterner = nullptr;
inline ASSEMBLYCPP_SEARCH_LOCAL std::unordered_map<
    treeCanonAtomId,
    treeCanonNodeId
> treeCanonLeafInternerDelta;
inline ASSEMBLYCPP_SEARCH_LOCAL std::uint64_t treeCanonInternerGeneration = 1;

void advanceTreeCanonInternerGeneration() noexcept
{
    if (++treeCanonInternerGeneration == 0) treeCanonInternerGeneration = 1;
}

void clearTreeCanonInterner()
{
    treeCanonInterner.clear();
    treeCanonAtomInterner.clear();
    treeCanonLeafInterner.clear();
    treeCanonLeafInternerDelta.clear();
    sharedTreeCanonAtomInterner = nullptr;
    sharedTreeCanonLeafInterner = nullptr;
    sharedTreeCanonInterner = nullptr;
    advanceTreeCanonInternerGeneration();
}

/** Borrow one immutable producer generation and retain only local misses. */
void bindTreeCanonInternerSeed(
    const decltype(treeCanonAtomInterner) &atomInterner,
    const decltype(treeCanonLeafInterner) &leafInterner,
    const decltype(treeCanonInterner) &treeInterner
)
{
    if (
        !treeCanonAtomInterner.empty() ||
        !treeCanonLeafInterner.empty() ||
        !treeCanonInterner.empty() ||
        !treeCanonLeafInternerDelta.empty()
    )
    {
        throw std::logic_error(
            "tree canonical seed requires empty worker deltas"
        );
    }
    sharedTreeCanonAtomInterner = &atomInterner;
    sharedTreeCanonLeafInterner = &leafInterner;
    sharedTreeCanonInterner = &treeInterner;
    advanceTreeCanonInternerGeneration();
}

treeCanonAtomId internTreeCanonAtom(const std::string &atomType)
{
    if (sharedTreeCanonAtomInterner == nullptr)
    {
        const treeCanonAtomId nextId =
            treeCanonAtomInterner.size() + 1;
        return treeCanonAtomInterner.try_emplace(
            atomType,
            nextId
        ).first->second;
    }
    const auto seeded = sharedTreeCanonAtomInterner->find(atomType);
    if (seeded != sharedTreeCanonAtomInterner->end()) return seeded->second;
    const treeCanonAtomId nextId =
        sharedTreeCanonAtomInterner->size() + treeCanonAtomInterner.size() + 1;
    return treeCanonAtomInterner.try_emplace(atomType, nextId).first->second;
}

treeCanonNodeId internTreeCanonSignature(treeCanonSignature signature)
{
    if (sharedTreeCanonInterner == nullptr)
    {
        const treeCanonNodeId nextId =
            treeCanonInterner.size() + 1;
        return treeCanonInterner.try_emplace(
            std::move(signature),
            nextId
        ).first->second;
    }
    // Seed IDs are dense and every worker-local ID starts above the seed
    // range. A signature containing any local atom or subtree therefore
    // cannot occur in the immutable seed; avoid a shared random lookup for
    // those overwhelmingly post-seed forms.
    bool canMatchSeed =
        signature.atomType <= sharedTreeCanonAtomInterner->size();
    if (canMatchSeed)
    {
        const treeCanonNodeId seededNodeCount =
            sharedTreeCanonInterner->size();
        for (const treeCanonChild &child : signature.children)
        {
            if (child.subtree > seededNodeCount)
            {
                canMatchSeed = false;
                break;
            }
        }
    }
    if (canMatchSeed)
    {
        const auto seeded = sharedTreeCanonInterner->find(signature);
        if (seeded != sharedTreeCanonInterner->end()) return seeded->second;
    }
    const treeCanonNodeId nextId =
        sharedTreeCanonInterner->size() + treeCanonInterner.size() + 1;
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
        if (sharedTreeCanonLeafInterner != nullptr)
        {
            if (
                atomType < sharedTreeCanonLeafInterner->size() &&
                (*sharedTreeCanonLeafInterner)[atomType] != 0
            )
            {
                return (*sharedTreeCanonLeafInterner)[atomType];
            }
            const auto localLeaf = treeCanonLeafInternerDelta.find(atomType);
            if (localLeaf != treeCanonLeafInternerDelta.end())
                return localLeaf->second;
            const treeCanonNodeId leaf =
                internTreeCanonSignature({atomType, {}});
            return treeCanonLeafInternerDelta.try_emplace(
                atomType,
                leaf
            ).first->second;
        }
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
    molGraph &graph,
    int node,
    std::vector<treeCanonChild> children
)
{
    return internTreeCanonNode(
        internTreeCanonAtom(graph.atoms[node].atomType),
        std::move(children)
    );
}

[[nodiscard]] inline std::size_t canonGraphVertexCount(
    const molGraph &graph
) noexcept
{
    return graph.atoms.size();
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
        graph.atoms.begin(),
        graph.atoms.end(),
        [](const atom &vertex) {return vertex.atomType == "X";}
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
        graph.atoms.begin(),
        graph.atoms.end(),
        [](const atom &vertex) {return vertex.bonds.size() <= 1;}
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
    return internTreeCanonAtom(graph.atoms[vertex].atomType);
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
    return graph.atoms[vertex].bonds;
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
    return edge.neighbourAtomIndex;
}

[[nodiscard]] inline int canonGraphNeighbour(
    const flatCanonAdjacentEdge &edge
) noexcept
{
    return static_cast<int>(edge.neighbour);
}

[[nodiscard]] inline std::uint16_t canonGraphBondType(const bond &edge) noexcept
{
    return static_cast<std::uint16_t>(edge.bondType);
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
 * @param graph Acyclic, connected molGraph
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

inline treeCanonForm centroidTreeCanon(molGraph &graph, int n)
{
    return centroidTreeCanonImpl(graph, n);
}

inline treeCanonForm centroidTreeCanon(const flatCanonGraph &graph, int n)
{
    return centroidTreeCanonImpl(graph, n);
}
