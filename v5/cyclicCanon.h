#pragma once

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

struct cyclicCanonVertexLabel
{
    // Keep raw atom colours and rooted-tree colours in distinct namespaces.
    std::uint64_t kind = 0;
    std::uint64_t value = 0;

    bool operator==(const cyclicCanonVertexLabel &) const = default;
};

using cyclicCanonAdjacentEdge = flatCanonAdjacentEdge;

struct cyclicCanonEdge
{
    std::uint32_t first = 0;
    std::uint32_t second = 0;
    std::uint16_t bondType = 0;
};

struct cyclicCanonGraph
{
    std::vector<cyclicCanonVertexLabel> labels;
    std::vector<std::size_t> adjacencyOffsets;
    std::vector<cyclicCanonAdjacentEdge> adjacency;
    std::vector<cyclicCanonEdge> edges;

    [[nodiscard]] std::span<const cyclicCanonAdjacentEdge> neighbours(
        std::size_t vertex
    ) const noexcept
    {
        return std::span<const cyclicCanonAdjacentEdge>(adjacency).subspan(
            adjacencyOffsets[vertex],
            adjacencyOffsets[vertex + 1] - adjacencyOffsets[vertex]
        );
    }
};

struct cyclicCanonPeelingWorkspace
{
    std::vector<int> degree;
    std::vector<int> localIndex;
    std::vector<unsigned char> removed;
    std::vector<std::vector<treeCanonChild>> children;
    std::vector<std::size_t> leaves;
    std::vector<std::size_t> coreVertices;
};

// The surrounding search state is already process-global and single-threaded.
// Reusing these flat buffers avoids repeated scratch allocation for every
// cyclic canonicalisation miss.
inline cyclicCanonPeelingWorkspace cyclicCanonPeelingScratch;

/**
 * @brief Cached coloured-core representation with lazy exact labelling.
 *
 * Construction peels pendant trees, materialises the core once, and stores a
 * deterministic integer-refinement hash used by unordered_map. Exact
 * individualisation/refinement is performed only when two candidates meet in
 * a hash bucket; its canonical code is then cached in the key.
 */
struct cyclicCanonForm
{
    cyclicCanonGraph graph;
    std::size_t cachedHash = 0;
    mutable std::vector<std::uint64_t> canonicalCode;

    [[nodiscard]] bool empty() const noexcept {return graph.labels.empty();}
    [[nodiscard]] std::size_t hash() const noexcept {return cachedHash;}
    [[nodiscard]] const std::vector<std::uint64_t> &exactCode() const;
    bool operator==(const cyclicCanonForm &other) const;
};

[[nodiscard]] bool cyclicCanonLabelLess(
    const cyclicCanonVertexLabel &left,
    const cyclicCanonVertexLabel &right
) noexcept
{
    if (left.kind != right.kind) return left.kind < right.kind;
    return left.value < right.value;
}

[[nodiscard]] std::uint16_t cyclicCanonBondType(short bondType) noexcept
{
    // Preserve the byte bond-label semantics used throughout canonicalisation.
    return static_cast<unsigned char>(static_cast<char>(bondType));
}

[[nodiscard]] std::uint32_t cyclicCanonColourCount(
    const std::vector<std::uint32_t> &colours
) noexcept
{
    if (colours.empty()) return 0;
    return *std::max_element(colours.begin(), colours.end()) + 1;
}

[[nodiscard]] std::vector<std::uint32_t> cyclicCanonInitialColours(
    const cyclicCanonGraph &graph
)
{
    const std::size_t vertexCount = graph.labels.size();
    if (vertexCount > std::numeric_limits<std::uint32_t>::max())
        throw std::length_error("cyclic canonical graph has too many vertices");

    std::vector<std::size_t> order(vertexCount);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(
        order.begin(),
        order.end(),
        [&](std::size_t left, std::size_t right)
        {
            return cyclicCanonLabelLess(graph.labels[left], graph.labels[right]);
        }
    );

    std::vector<std::uint32_t> colours(vertexCount, 0);
    std::uint32_t colour = 0;
    for (std::size_t position = 0; position < order.size(); position++)
    {
        if (
            position != 0 &&
            !(graph.labels[order[position]] == graph.labels[order[position - 1]])
        ) ++colour;
        colours[order[position]] = colour;
    }
    return colours;
}

/**
 * @brief Deterministic 1-dimensional Weisfeiler-Leman colour refinement.
 *
 * Integer colours are assigned by sorting exact signatures, never by hashing,
 * so refinement cannot merge different neighbourhood multisets by collision.
 */
[[nodiscard]] std::vector<std::uint32_t> cyclicCanonRefineColours(
    const cyclicCanonGraph &graph,
    std::vector<std::uint32_t> colours
)
{
    const std::size_t vertexCount = graph.labels.size();
    std::uint32_t previousCount = cyclicCanonColourCount(colours);
    if (previousCount == vertexCount) return colours;

    // Store every sorted (bond, neighbour-colour) token in one buffer. This
    // retains exact multiset comparison without one heap allocation per core
    // vertex and per refinement call.
    const std::vector<std::size_t> &offsets = graph.adjacencyOffsets;
    std::vector<std::uint64_t> neighbourTokens(graph.adjacency.size());
    std::vector<std::size_t> order(vertexCount);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::vector<std::uint32_t> refined(vertexCount, 0);

    while (true)
    {
        for (std::size_t vertex = 0; vertex < vertexCount; vertex++)
        {
            std::size_t token = offsets[vertex];
            for (const cyclicCanonAdjacentEdge &edge : graph.neighbours(vertex))
            {
                neighbourTokens[token++] =
                    (static_cast<std::uint64_t>(edge.bondType) << 32) |
                    colours[edge.neighbour];
            }
            std::sort(
                neighbourTokens.begin() + offsets[vertex],
                neighbourTokens.begin() + offsets[vertex + 1]
            );
        }

        std::sort(
            order.begin(),
            order.end(),
            [&](std::size_t left, std::size_t right)
            {
                if (colours[left] != colours[right])
                    return colours[left] < colours[right];
                return std::lexicographical_compare(
                    neighbourTokens.begin() + offsets[left],
                    neighbourTokens.begin() + offsets[left + 1],
                    neighbourTokens.begin() + offsets[right],
                    neighbourTokens.begin() + offsets[right + 1]
                );
            }
        );

        std::uint32_t refinedColour = 0;
        for (std::size_t position = 0; position < order.size(); position++)
        {
            if (position != 0)
            {
                const std::size_t current = order[position];
                const std::size_t previous = order[position - 1];
                const bool equal =
                    colours[current] == colours[previous] &&
                    std::equal(
                        neighbourTokens.begin() + offsets[current],
                        neighbourTokens.begin() + offsets[current + 1],
                        neighbourTokens.begin() + offsets[previous],
                        neighbourTokens.begin() + offsets[previous + 1]
                    );
                if (!equal) ++refinedColour;
            }
            refined[order[position]] = refinedColour;
        }
        const std::uint32_t refinedCount = vertexCount == 0
            ? 0 : refinedColour + 1;
        if (refinedCount == previousCount || refinedCount == vertexCount)
            return refined;
        colours.swap(refined);
        previousCount = refinedCount;
    }
}

[[nodiscard]] std::vector<std::uint64_t> cyclicCanonEncode(
    const cyclicCanonGraph &graph,
    const std::vector<std::uint32_t> &discreteColours
)
{
    const std::size_t vertexCount = graph.labels.size();
    std::vector<std::uint32_t> vertexAtColour(vertexCount, 0);
    for (std::size_t vertex = 0; vertex < vertexCount; vertex++)
        vertexAtColour[discreteColours[vertex]] = static_cast<std::uint32_t>(vertex);

    std::vector<std::uint64_t> code;
    code.reserve(3 + vertexCount * 2 + graph.edges.size() * 3);
    code.push_back(0); // General coloured-core canonical code.
    code.push_back(vertexCount);
    code.push_back(graph.edges.size());
    for (const std::uint32_t vertex : vertexAtColour)
    {
        code.push_back(graph.labels[vertex].kind);
        code.push_back(graph.labels[vertex].value);
    }

    std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint16_t>> edges;
    edges.reserve(graph.edges.size());
    for (const cyclicCanonEdge &edge : graph.edges)
    {
        std::uint32_t first = discreteColours[edge.first];
        std::uint32_t second = discreteColours[edge.second];
        if (second < first) std::swap(first, second);
        edges.emplace_back(first, second, edge.bondType);
    }
    std::sort(edges.begin(), edges.end());
    for (const auto &[first, second, bondType] : edges)
    {
        code.push_back(first);
        code.push_back(second);
        code.push_back(bondType);
    }
    return code;
}

[[nodiscard]] bool cyclicCanonIsSimpleCycle(
    const cyclicCanonGraph &graph
)
{
    if (graph.labels.size() < 3 || graph.edges.size() != graph.labels.size())
        return false;
    for (std::size_t vertex = 0; vertex < graph.labels.size(); vertex++)
    {
        const auto adjacency = graph.neighbours(vertex);
        if (
            adjacency.size() != 2 ||
            adjacency[0].neighbour == adjacency[1].neighbour
        ) return false;
    }

    // Degree two and |E| == |V| also describes a disjoint union of cycles.
    // The dihedral encoder is valid only when the whole representation is one
    // connected cycle.
    std::vector<unsigned char> seen(graph.labels.size(), 0);
    std::vector<std::uint32_t> pending{0};
    seen[0] = 1;
    std::size_t reached = 0;
    while (!pending.empty())
    {
        const std::uint32_t vertex = pending.back();
        pending.pop_back();
        ++reached;
        for (const cyclicCanonAdjacentEdge &edge : graph.neighbours(vertex))
        {
            if (edge.neighbour >= graph.labels.size()) return false;
            if (seen[edge.neighbour]) continue;
            seen[edge.neighbour] = 1;
            pending.push_back(edge.neighbour);
        }
    }
    return reached == graph.labels.size();
}

/** Exact dihedral canonical form for the common one-cycle 2-core. */
[[nodiscard]] std::vector<std::uint64_t> cyclicCanonEncodeSimpleCycle(
    const cyclicCanonGraph &graph
)
{
    const std::size_t vertexCount = graph.labels.size();
    std::vector<std::uint64_t> best;
    for (std::size_t start = 0; start < vertexCount; start++)
    {
        for (std::size_t direction = 0; direction < 2; direction++)
        {
            std::vector<std::uint64_t> candidate;
            candidate.reserve(2 + vertexCount * 3);
            candidate.push_back(1); // Simple-cycle canonical code.
            candidate.push_back(vertexCount);

            std::size_t previous = vertexCount;
            std::size_t current = start;
            std::size_t next = graph.neighbours(start)[direction].neighbour;
            for (std::size_t step = 0; step < vertexCount; step++)
            {
                candidate.push_back(graph.labels[current].kind);
                candidate.push_back(graph.labels[current].value);

                std::uint16_t bondType = 0;
                bool bondFound = false;
                for (const cyclicCanonAdjacentEdge &edge : graph.neighbours(current))
                {
                    if (edge.neighbour == next)
                    {
                        bondType = edge.bondType;
                        bondFound = true;
                        break;
                    }
                }
                if (!bondFound)
                    throw std::logic_error("cyclic canonical cycle edge is invalid");
                candidate.push_back(bondType);

                previous = current;
                current = next;
                if (step + 1 < vertexCount)
                {
                    const auto adjacency = graph.neighbours(current);
                    next = adjacency[0].neighbour == previous
                        ? adjacency[1].neighbour
                        : adjacency[0].neighbour;
                }
            }
            if (current != start)
                throw std::logic_error("cyclic canonical core is not a cycle");
            if (best.empty() || candidate < best) best = std::move(candidate);
        }
    }
    return best;
}

class cyclicCanonSearch
{
public:
    explicit cyclicCanonSearch(const cyclicCanonGraph &_graph): graph(_graph) {}

    [[nodiscard]] std::vector<std::uint64_t> run()
    {
        visit(cyclicCanonInitialColours(graph));
        return std::move(bestCode);
    }

    [[nodiscard]] std::vector<std::uint64_t> run(
        std::vector<std::uint32_t> initialColours
    )
    {
        visit(std::move(initialColours));
        return std::move(bestCode);
    }

private:
    void visit(std::vector<std::uint32_t> colours)
    {
        colours = cyclicCanonRefineColours(graph, std::move(colours));
        const std::uint32_t colourCount = cyclicCanonColourCount(colours);
        if (colourCount == graph.labels.size())
        {
            std::vector<std::uint64_t> code = cyclicCanonEncode(graph, colours);
            if (bestCode.empty() || code < bestCode) bestCode = std::move(code);
            return;
        }

        std::vector<std::size_t> cellSizes(colourCount, 0);
        for (const std::uint32_t colour : colours) ++cellSizes[colour];
        std::size_t selectedSize = std::numeric_limits<std::size_t>::max();
        std::uint32_t selectedColour = 0;
        for (std::uint32_t colour = 0; colour < colourCount; colour++)
        {
            if (cellSizes[colour] > 1 && cellSizes[colour] < selectedSize)
            {
                selectedSize = cellSizes[colour];
                selectedColour = colour;
            }
        }

        for (std::size_t vertex = 0; vertex < colours.size(); vertex++)
        {
            if (colours[vertex] != selectedColour) continue;
            std::vector<std::uint32_t> individualized = colours;
            individualized[vertex] = colourCount;
            visit(std::move(individualized));
        }
    }

    const cyclicCanonGraph &graph;
    std::vector<std::uint64_t> bestCode;
};

[[nodiscard]] std::size_t cyclicCanonHashCode(
    const std::vector<std::uint64_t> &code
) noexcept
{
    std::size_t result = static_cast<std::size_t>(0xcbf29ce484222325ULL);
    for (const std::uint64_t value : code)
    {
        treeCanonSignatureHash::combine(
            result,
            std::hash<std::uint64_t>{}(value)
        );
    }
    return result;
}

[[nodiscard]] std::uint64_t cyclicCanonMix(std::uint64_t value) noexcept
{
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

[[nodiscard]] std::size_t cyclicCanonCoarseHash(
    const cyclicCanonGraph &graph,
    const std::vector<std::uint32_t> &stableColours
) noexcept
{
    const std::uint32_t colourCount = cyclicCanonColourCount(stableColours);
    std::uint64_t vertexSum = 0;
    std::uint64_t vertexXor = 0;
    for (std::size_t vertex = 0; vertex < graph.labels.size(); vertex++)
    {
        std::uint64_t vertexHash = cyclicCanonMix(stableColours[vertex]);
        vertexHash = cyclicCanonMix(vertexHash ^ graph.labels[vertex].kind);
        vertexHash = cyclicCanonMix(vertexHash ^ graph.labels[vertex].value);
        vertexSum += vertexHash;
        vertexXor ^= std::rotl(
            vertexHash,
            static_cast<int>(vertexHash & 63)
        );
    }

    std::uint64_t result = cyclicCanonMix(graph.labels.size());
    result = cyclicCanonMix(result ^ graph.edges.size());
    result = cyclicCanonMix(result ^ colourCount);
    result = cyclicCanonMix(result ^ vertexSum);
    result = cyclicCanonMix(result ^ vertexXor);

    std::uint64_t edgeSum = 0;
    std::uint64_t edgeXor = 0;
    for (const cyclicCanonEdge &edge : graph.edges)
    {
        std::uint32_t first = stableColours[edge.first];
        std::uint32_t second = stableColours[edge.second];
        if (second < first) std::swap(first, second);
        std::uint64_t edgeHash = cyclicCanonMix(first);
        edgeHash = cyclicCanonMix(edgeHash ^ second);
        edgeHash = cyclicCanonMix(edgeHash ^ edge.bondType);
        edgeSum += edgeHash;
        edgeXor ^= std::rotl(edgeHash, static_cast<int>(edgeHash & 63));
    }
    result = cyclicCanonMix(result ^ edgeSum);
    result = cyclicCanonMix(result ^ edgeXor);
    return static_cast<std::size_t>(result);
}

const std::vector<std::uint64_t> &cyclicCanonForm::exactCode() const
{
    if (!canonicalCode.empty()) return canonicalCode;
    if (cyclicCanonIsSimpleCycle(graph))
    {
        canonicalCode = cyclicCanonEncodeSimpleCycle(graph);
    }
    else
    {
        cyclicCanonSearch search(graph);
        canonicalCode = search.run();
    }
    return canonicalCode;
}

bool cyclicCanonForm::operator==(const cyclicCanonForm &other) const
{
    return exactCode() == other.exactCode();
}

[[nodiscard]] cyclicCanonForm canonicaliseCyclicCanonGraph(
    cyclicCanonGraph graph
)
{
    cyclicCanonForm result;
    std::vector<std::uint32_t> stableColours = cyclicCanonRefineColours(
        graph,
        cyclicCanonInitialColours(graph)
    );
    result.cachedHash = cyclicCanonCoarseHash(graph, stableColours);
    result.graph = std::move(graph);
    return result;
}

template<typename Graph>
[[nodiscard]] cyclicCanonGraph buildWholeGraphCanonRepresentationImpl(
    const Graph &input,
    std::uint64_t labelKind = 0
)
{
    const std::size_t vertexCount = canonGraphVertexCount(input);
    if (vertexCount > std::numeric_limits<std::uint32_t>::max())
        throw std::length_error("cyclic canonical graph has too many vertices");

    cyclicCanonGraph graph;
    graph.labels.reserve(vertexCount);
    graph.adjacencyOffsets.reserve(vertexCount + 1);
    graph.adjacencyOffsets.push_back(0);
    std::size_t adjacencyCount = 0;
    for (std::size_t vertex = 0; vertex < vertexCount; vertex++)
    {
        graph.labels.push_back({
            labelKind,
            canonGraphAtomType(input, vertex)
        });
        adjacencyCount += canonGraphNeighbours(input, vertex).size();
    }
    graph.adjacency.reserve(adjacencyCount);
    graph.edges.reserve(adjacencyCount / 2);
    for (std::size_t first = 0; first < vertexCount; first++)
    {
        for (const auto &edge : canonGraphNeighbours(input, first))
        {
            const int neighbour = canonGraphNeighbour(edge);
            if (neighbour < 0 || static_cast<std::size_t>(neighbour) >= vertexCount)
                throw std::logic_error("cyclic canonical graph edge is invalid");
            const std::size_t second = static_cast<std::size_t>(neighbour);
            const std::uint16_t bondType = canonGraphBondType(edge);
            graph.adjacency.push_back({
                static_cast<std::uint32_t>(second), bondType
            });
            if (first < second)
            {
                graph.edges.push_back({
                    static_cast<std::uint32_t>(first),
                    static_cast<std::uint32_t>(second),
                    bondType
                });
            }
        }
        graph.adjacencyOffsets.push_back(graph.adjacency.size());
    }
    return graph;
}

/** Bulk-copy an already-flat input into the owned lazy canonical form. */
[[nodiscard]] cyclicCanonGraph buildWholeGraphCanonRepresentationImpl(
    const flatCanonGraph &input,
    std::uint64_t labelKind = 0
)
{
    const std::size_t vertexCount = input.labels.size();
    if (vertexCount > std::numeric_limits<std::uint32_t>::max())
        throw std::length_error("cyclic canonical graph has too many vertices");

    cyclicCanonGraph graph;
    graph.labels.resize(vertexCount);
    for (std::size_t vertex = 0; vertex < vertexCount; vertex++)
        graph.labels[vertex] = {labelKind, input.labels[vertex]};
    graph.adjacencyOffsets = input.adjacencyOffsets;
    graph.adjacency = input.adjacency;
    graph.edges.reserve(input.edgeCount);
    for (std::size_t first = 0; first < vertexCount; first++)
    {
        for (const flatCanonAdjacentEdge &edge : input.neighbours(first))
        {
            if (first < edge.neighbour)
            {
                graph.edges.push_back({
                    static_cast<std::uint32_t>(first),
                    edge.neighbour,
                    edge.bondType
                });
            }
        }
    }
    return graph;
}

/**
 * @brief Build the exact coloured 2-core of a connected cyclic graph.
 *
 * Pendant trees are peeled iteratively and interned using the same exact AHU
 * signatures as the acyclic canonicaliser. Their rooted identity becomes the
 * integer colour of the core vertex to which they are attached.
 */
template<typename Graph>
[[nodiscard]] cyclicCanonGraph buildColouredTwoCoreRepresentationImpl(
    const Graph &input
)
{
    const std::size_t vertexCount = canonGraphVertexCount(input);
    if (vertexCount > std::numeric_limits<std::uint32_t>::max())
        throw std::length_error("cyclic canonical graph has too many vertices");

    // The common ring/bridged-ring case is already its own 2-core. Avoid all
    // peeling scratch and rooted-tree bookkeeping when no pendant vertex can
    // be removed.
    const bool requiresPeeling = canonGraphHasPendantVertex(input);
    if (!requiresPeeling)
        return buildWholeGraphCanonRepresentationImpl(input, 1);

    auto &degree = cyclicCanonPeelingScratch.degree;
    auto &removed = cyclicCanonPeelingScratch.removed;
    auto &children = cyclicCanonPeelingScratch.children;
    auto &leaves = cyclicCanonPeelingScratch.leaves;
    degree.resize(vertexCount);
    removed.assign(vertexCount, 0);
    children.resize(vertexCount);
    for (auto &vertexChildren : children) vertexChildren.clear();
    leaves.clear();
    if (leaves.capacity() < vertexCount) leaves.reserve(vertexCount);
    for (std::size_t vertex = 0; vertex < vertexCount; vertex++)
    {
        degree[vertex] = static_cast<int>(
            canonGraphNeighbours(input, vertex).size()
        );
        if (degree[vertex] <= 1) leaves.push_back(vertex);
    }

    bool detachedTreeComponent = false;
    for (std::size_t position = 0; position < leaves.size(); position++)
    {
        const std::size_t leaf = leaves[position];
        if (removed[leaf]) continue;
        removed[leaf] = 1;
        const treeCanonNodeId subtree = internTreeCanonNode(
            canonGraphAtomType(input, leaf),
            std::move(children[leaf])
        );

        std::size_t liveNeighbour = vertexCount;
        std::uint16_t liveBond = 0;
        for (const auto &edge : canonGraphNeighbours(input, leaf))
        {
            const int neighbourIndex = canonGraphNeighbour(edge);
            if (
                neighbourIndex < 0 ||
                static_cast<std::size_t>(neighbourIndex) >= vertexCount
            )
                throw std::logic_error("cyclic canonical graph edge is invalid");
            const std::size_t neighbour = static_cast<std::size_t>(neighbourIndex);
            if (removed[neighbour]) continue;
            liveNeighbour = neighbour;
            liveBond = canonGraphBondType(edge);
            break;
        }
        if (liveNeighbour == vertexCount)
        {
            detachedTreeComponent = true;
            continue;
        }
        children[liveNeighbour].push_back({
            static_cast<unsigned char>(liveBond),
            subtree
        });
        if (--degree[liveNeighbour] == 1) leaves.push_back(liveNeighbour);
    }

    auto &localIndex = cyclicCanonPeelingScratch.localIndex;
    auto &coreVertices = cyclicCanonPeelingScratch.coreVertices;
    localIndex.assign(vertexCount, -1);
    coreVertices.clear();
    if (coreVertices.capacity() < vertexCount) coreVertices.reserve(vertexCount);
    std::size_t coreAdjacencyCount = 0;
    for (std::size_t vertex = 0; vertex < vertexCount; vertex++)
    {
        if (removed[vertex]) continue;
        localIndex[vertex] = static_cast<int>(coreVertices.size());
        coreVertices.push_back(vertex);
        coreAdjacencyCount += static_cast<std::size_t>(degree[vertex]);
    }
    if (coreVertices.empty() || detachedTreeComponent)
        return buildWholeGraphCanonRepresentationImpl(input);

    cyclicCanonGraph graph;
    graph.labels.reserve(coreVertices.size());
    graph.adjacencyOffsets.reserve(coreVertices.size() + 1);
    graph.adjacencyOffsets.push_back(0);
    graph.adjacency.reserve(coreAdjacencyCount);
    graph.edges.reserve(coreAdjacencyCount / 2);
    for (const std::size_t vertex : coreVertices)
    {
        if (children[vertex].empty())
        {
            graph.labels.push_back({
                1,
                canonGraphAtomType(input, vertex)
            });
        }
        else
        {
            graph.labels.push_back({
                2,
                internTreeCanonNode(
                    canonGraphAtomType(input, vertex),
                    std::move(children[vertex])
                )
            });
        }
    }

    for (std::size_t localFirst = 0; localFirst < coreVertices.size(); localFirst++)
    {
        const std::size_t first = coreVertices[localFirst];
        for (const auto &edge : canonGraphNeighbours(input, first))
        {
            const int neighbour = canonGraphNeighbour(edge);
            if (neighbour < 0 || static_cast<std::size_t>(neighbour) >= vertexCount)
                throw std::logic_error("cyclic canonical core edge is invalid");
            const std::size_t second = static_cast<std::size_t>(neighbour);
            if (removed[second]) continue;
            const int localSecond = localIndex[second];
            if (localSecond < 0)
                throw std::logic_error("cyclic canonical core edge is invalid");
            const std::uint16_t bondType = canonGraphBondType(edge);
            graph.adjacency.push_back({
                static_cast<std::uint32_t>(localSecond), bondType
            });
            if (localFirst < static_cast<std::size_t>(localSecond))
            {
                graph.edges.push_back({
                    static_cast<std::uint32_t>(localFirst),
                    static_cast<std::uint32_t>(localSecond),
                    bondType
                });
            }
        }
        graph.adjacencyOffsets.push_back(graph.adjacency.size());
    }
    return graph;
}

[[nodiscard]] inline cyclicCanonGraph buildWholeGraphCanonRepresentation(
    molGraph &mg,
    std::uint64_t labelKind = 0
)
{
    return buildWholeGraphCanonRepresentationImpl(mg, labelKind);
}

[[nodiscard]] inline cyclicCanonGraph buildWholeGraphCanonRepresentation(
    const flatCanonGraph &graph,
    std::uint64_t labelKind = 0
)
{
    return buildWholeGraphCanonRepresentationImpl(graph, labelKind);
}

[[nodiscard]] inline cyclicCanonGraph buildColouredTwoCoreRepresentation(
    molGraph &mg
)
{
    return buildColouredTwoCoreRepresentationImpl(mg);
}

[[nodiscard]] inline cyclicCanonGraph buildColouredTwoCoreRepresentation(
    const flatCanonGraph &graph
)
{
    return buildColouredTwoCoreRepresentationImpl(graph);
}

[[nodiscard]] inline cyclicCanonForm canonicaliseCyclicGraph(molGraph &mg)
{
    return canonicaliseCyclicCanonGraph(buildColouredTwoCoreRepresentation(mg));
}

[[nodiscard]] inline cyclicCanonForm canonicaliseCyclicGraph(
    const flatCanonGraph &graph
)
{
    return canonicaliseCyclicCanonGraph(
        buildColouredTwoCoreRepresentation(graph)
    );
}

[[nodiscard]] inline cyclicCanonForm canonicaliseWholeGraph(molGraph &mg)
{
    return canonicaliseCyclicCanonGraph(buildWholeGraphCanonRepresentation(mg));
}

[[nodiscard]] inline cyclicCanonForm canonicaliseWholeGraph(
    const flatCanonGraph &graph
)
{
    return canonicaliseCyclicCanonGraph(
        buildWholeGraphCanonRepresentation(graph)
    );
}
