#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <mutex>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include "compilerAttributes.h"

/**
 * @brief Hashes a molecular graph
 */
struct graphHash
{
    /// Exact interned tree identity; empty selects the exact general form.
    treeCanonForm treeHash;
    /// Cached coloured-core or whole-graph form with a lazy exact code.
    cyclicCanonForm cyclicHash;

    /**
     * @brief Construct a new graph Hash object
     *
     * @param graph molGraph to be hashed
     * @param isCyclic Is the molecule cyclic
     */
    graphHash(molGraph &graph, bool isCyclic)
    {
        if (isCyclic) cyclicHash = canonicaliseCyclicGraph(graph);
        else
        {
            treeHash = centroidTreeCanon(graph, 0);
            if (treeHash.empty()) cyclicHash = canonicaliseWholeGraph(graph);
        }
    }

    /** Build a self-contained hash directly from reusable flat CSR input. */
    graphHash(const flatCanonGraph &graph, bool isCyclic)
    {
        if (isCyclic) cyclicHash = canonicaliseCyclicGraph(graph);
        else
        {
            treeHash = centroidTreeCanon(graph, 0);
            if (treeHash.empty()) cyclicHash = canonicaliseWholeGraph(graph);
        }
    }

    /**
     * @brief Compare cached exact canonical forms.
     *
     * @param g2 other graph to be compared
     * @return true if graphs are isomorphic
     * @return false otherwise
     */
    bool operator==(const graphHash &g2) const
    {
        if (treeHash.empty() != g2.treeHash.empty()) return false;
        if (treeHash.empty()) return cyclicHash == g2.cyclicHash;
        return treeHash == g2.treeHash;
    }

    /** Materialise the only lazily mutable field before sharing this key. */
    void prepareForSharing() const
    {
        if (treeHash.empty()) static_cast<void>(cyclicHash.exactCode());
    }
};

/**
 * @brief Hash for subgraph unordered_map
 */
#ifdef ASSEMBLYCPP_LIBRARY_BUILD
struct graphHashHasher
#else
template<>
struct std::hash<graphHash>
#endif
{
    size_t operator()(const graphHash &gh) const
    {
        if (gh.treeHash.empty())
        {
            return gh.cyclicHash.hash();
        }
        else
        {
            size_t result = std::hash<treeCanonNodeId>{}(gh.treeHash.first);
            treeCanonSignatureHash::combine(
                result,
                std::hash<treeCanonNodeId>{}(gh.treeHash.second)
            );
            treeCanonSignatureHash::combine(
                result,
                std::hash<std::uint16_t>{}(gh.treeHash.centralBond)
            );
            return result;
        }
    }
};

/**
 * @brief Reusable edge-mask to canonical CSR builder.
 *
 * Source atom and bond labels are interned/normalized once per calculation.
 * Each miss then visits only selected mask bits, maps their incident vertices,
 * and fills flat adjacency without constructing atoms, strings, or one vector
 * per vertex.
 */
struct canonicalisationGraphWorkspace
{
    struct flatEdge
    {
        std::uint32_t first = 0;
        std::uint32_t second = 0;
        std::uint16_t bondType = 0;
    };

    std::vector<treeCanonAtomId> sourceLabels;
    std::vector<unsigned char> sourceLegacyX;
    std::vector<flatEdge> sourceEdges;
    std::vector<int> localVertex;
    std::vector<std::uint32_t> vertexEpoch;
    std::uint32_t epoch = 0;
    std::vector<disjointSetNode> unionNodes;
    std::vector<std::size_t> degree;
    std::vector<flatEdge> selectedEdges;
    flatCanonGraph graph;
    std::uint64_t configuredInternerGeneration = 0;
    const atom *configuredAtomData = nullptr;
    const MoleculeEdge *configuredEdgeData = nullptr;

    /** Must be called after changing source atoms or the universe edge list. */
    void configure(
        const molGraph &source,
        const std::vector<MoleculeEdge> &edgeList
    )
    {
        if (source.atoms.size() > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("canonical graph has too many vertices");

        sourceLabels.clear();
        sourceLabels.reserve(source.atoms.size());
        sourceLegacyX.clear();
        sourceLegacyX.reserve(source.atoms.size());
        for (const atom &vertex : source.atoms)
        {
            sourceLabels.push_back(internTreeCanonAtom(vertex.atomType));
            sourceLegacyX.push_back(vertex.atomType == "X");
        }

        sourceEdges.clear();
        sourceEdges.reserve(edgeList.size());
        for (const MoleculeEdge &edge : edgeList)
        {
            if (
                edge.sourceAtomIndex < 0 || edge.targetAtomIndex < 0 ||
                edge.sourceBondIndex < 0 ||
                edge.sourceAtomIndex == edge.targetAtomIndex ||
                static_cast<std::size_t>(edge.sourceAtomIndex) >=
                    source.atoms.size() ||
                static_cast<std::size_t>(edge.targetAtomIndex) >=
                    source.atoms.size() ||
                static_cast<std::size_t>(edge.sourceBondIndex) >=
                    source.atoms[
                        static_cast<std::size_t>(edge.sourceAtomIndex)
                    ].bonds.size()
            )
            {
                throw std::logic_error("canonical source edge is invalid");
            }
            const bond &sourceBond =
                source.atoms[edge.sourceAtomIndex]
                    .bonds[edge.sourceBondIndex];
            sourceEdges.push_back({
                static_cast<std::uint32_t>(edge.sourceAtomIndex),
                static_cast<std::uint32_t>(edge.targetAtomIndex),
                canonGraphBondType(sourceBond)
            });
        }

        localVertex.resize(source.atoms.size());
        vertexEpoch.resize(source.atoms.size(), 0);
        unionNodes.reserve(source.atoms.size());
        degree.reserve(source.atoms.size());
        selectedEdges.reserve(edgeList.size());
        graph.labels.reserve(source.atoms.size());
        graph.adjacencyOffsets.reserve(source.atoms.size() + 1);
        graph.adjacency.reserve(edgeList.size() * 2);
        configuredInternerGeneration = treeCanonInternerGeneration;
        configuredAtomData = source.atoms.data();
        configuredEdgeData = edgeList.data();
    }

    void ensureConfigured(
        const molGraph &source,
        const std::vector<MoleculeEdge> &edgeList
    )
    {
        // This catches interner resets and storage replacement. Call configure
        // explicitly after any in-place atom or edge content mutation, even
        // when vector sizes and backing addresses are unchanged.
        if (
            configuredInternerGeneration != treeCanonInternerGeneration ||
            sourceLabels.size() != source.atoms.size() ||
            sourceEdges.size() != edgeList.size() ||
            configuredAtomData != source.atoms.data() ||
            configuredEdgeData != edgeList.data()
        )
        {
            configure(source, edgeList);
        }
    }

    void begin()
    {
        if (++epoch == 0)
        {
            std::fill(vertexEpoch.begin(), vertexEpoch.end(), 0);
            epoch = 1;
        }
        unionNodes.clear();
        degree.clear();
        selectedEdges.clear();
        graph.labels.clear();
        graph.adjacencyOffsets.clear();
        graph.adjacency.clear();
        graph.edgeCount = 0;
        graph.hasLegacyX = false;
        graph.hasPendantVertex = false;
    }

    [[nodiscard]] bool contains(std::size_t sourceVertex) const noexcept
    {
        return vertexEpoch[sourceVertex] == epoch;
    }

    int addVertex(std::size_t sourceVertex, int parent = -1)
    {
        const std::size_t local = unionNodes.size();
        if (local > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("canonical graph has too many vertices");
        localVertex[sourceVertex] = static_cast<int>(local);
        vertexEpoch[sourceVertex] = epoch;
        unionNodes.emplace_back();
        unionNodes.back().parent = parent < 0
            ? static_cast<int>(local) : parent;
        unionNodes.back().rank = 0;
        degree.push_back(0);
        graph.labels.push_back(sourceLabels[sourceVertex]);
        graph.hasLegacyX |= sourceLegacyX[sourceVertex] != 0;
        return static_cast<int>(local);
    }

    [[nodiscard]] std::size_t find(std::size_t vertex)
    {
        disjointSetNode &node = unionNodes[vertex];
        if (node.parent != static_cast<int>(vertex))
            node.parent = static_cast<int>(find(node.parent));
        return static_cast<std::size_t>(node.parent);
    }

    /** Return true when the edge closes a cycle. */
    bool merge(std::size_t first, std::size_t second)
    {
        const std::size_t firstRoot = find(first);
        const std::size_t secondRoot = find(second);
        if (firstRoot == secondRoot) return true;
        if (unionNodes[firstRoot].rank > unionNodes[secondRoot].rank)
        {
            unionNodes[secondRoot].parent = static_cast<int>(firstRoot);
        }
        else
        {
            unionNodes[firstRoot].parent = static_cast<int>(secondRoot);
            if (unionNodes[firstRoot].rank == unionNodes[secondRoot].rank)
                ++unionNodes[secondRoot].rank;
        }
        return false;
    }

    /** Add one selected source edge without duplicating this large body in
     * the one-word and wide-mask iteration paths. */
    ASSEMBLYCPP_NOINLINE void addSelectedEdge(
        std::size_t edgeIndex,
        bool &isCyclic
    )
    {
        const flatEdge &sourceEdge = sourceEdges[edgeIndex];
        const std::size_t sourceFirst = sourceEdge.first;
        const std::size_t sourceSecond = sourceEdge.second;
        const bool hasFirst = contains(sourceFirst);
        const bool hasSecond = contains(sourceSecond);
        int first;
        int second;
        if (!hasFirst && !hasSecond)
        {
            first = addVertex(sourceFirst);
            second = addVertex(sourceSecond, first);
        }
        else if (!hasFirst)
        {
            second = localVertex[sourceSecond];
            first = addVertex(sourceFirst, second);
        }
        else if (!hasSecond)
        {
            first = localVertex[sourceFirst];
            second = addVertex(sourceSecond, first);
        }
        else
        {
            first = localVertex[sourceFirst];
            second = localVertex[sourceSecond];
            if (!isCyclic)
            {
                isCyclic = merge(
                    static_cast<std::size_t>(first),
                    static_cast<std::size_t>(second)
                );
            }
        }
        ++degree[static_cast<std::size_t>(first)];
        ++degree[static_cast<std::size_t>(second)];
        selectedEdges.push_back({
            static_cast<std::uint32_t>(first),
            static_cast<std::uint32_t>(second),
            sourceEdge.bondType
        });
    }

    const flatCanonGraph &build(
        const molGraph &source,
        const std::vector<MoleculeEdge> &edgeList,
        const EdgeMask &mask,
        bool &isCyclic
    )
    {
        ensureConfigured(source, edgeList);
        begin();
        isCyclic = false;
        forEachSetBitBelow(mask, sourceEdges.size(), [&](std::size_t edgeIndex)
        {
            addSelectedEdge(edgeIndex, isCyclic);
        });

        graph.edgeCount = selectedEdges.size();
        graph.adjacencyOffsets.resize(graph.labels.size() + 1);
        graph.adjacencyOffsets[0] = 0;
        for (std::size_t vertex = 0; vertex < graph.labels.size(); vertex++)
        {
            graph.hasPendantVertex |= degree[vertex] <= 1;
            graph.adjacencyOffsets[vertex + 1] =
                graph.adjacencyOffsets[vertex] + degree[vertex];
            degree[vertex] = graph.adjacencyOffsets[vertex];
        }
        graph.adjacency.resize(selectedEdges.size() * 2);
        for (const flatEdge &edge : selectedEdges)
        {
            graph.adjacency[degree[edge.first]++] = {
                edge.second,
                edge.bondType
            };
            graph.adjacency[degree[edge.second]++] = {
                edge.first,
                edge.bondType
            };
        }
        return graph;
    }
};

/** A cheap, deterministic isomorphism invariant used only to select a shard
 * bucket. Exact graphHash comparison still decides equality. */
[[nodiscard]] std::size_t sharedCanonicalInvariantHash(
    const flatCanonGraph &graph
) noexcept
{
    auto mix = [](std::uint64_t value) noexcept
    {
        value ^= value >> 30;
        value *= UINT64_C(0xbf58476d1ce4e5b9);
        value ^= value >> 27;
        value *= UINT64_C(0x94d049bb133111eb);
        value ^= value >> 31;
        return value;
    };

    std::uint64_t vertexSum = 0;
    std::uint64_t vertexXor = 0;
    for (std::size_t vertex = 0; vertex < graph.labels.size(); ++vertex)
    {
        std::uint64_t bondSum = 0;
        std::uint64_t bondXor = 0;
        const auto neighbours = graph.neighbours(vertex);
        for (const flatCanonAdjacentEdge &edge : neighbours)
        {
            const std::uint64_t token = mix(edge.bondType);
            bondSum += token;
            bondXor ^= std::rotl(
                token,
                static_cast<int>(token & UINT64_C(63))
            );
        }
        std::uint64_t token = mix(graph.labels[vertex]);
        token = mix(token ^ neighbours.size());
        token = mix(token ^ bondSum);
        token = mix(token ^ bondXor);
        vertexSum += token;
        vertexXor ^= std::rotl(
            token,
            static_cast<int>(token & UINT64_C(63))
        );
    }

    std::uint64_t result = mix(graph.labels.size());
    result = mix(result ^ graph.edgeCount);
    result = mix(result ^ vertexSum);
    result = mix(result ^ vertexXor);
    return static_cast<std::size_t>(result);
}

/**
 * Process-shared canonical-ID allocator with exact out-of-lock comparisons.
 *
 * The shared representatives are serialized source masks, not worker-owned
 * graphHash objects: tree and peeled-core forms contain worker-local interner
 * IDs. A worker snapshots one immutable representative under a shard lock,
 * reconstructs and canonicalises it in its own interner generation after
 * releasing the lock, and returns only after exact graphHash equality. If no
 * representative matches, insertion is checked again under the lock. Every
 * worker first checks the same immutable producer graph-hash seed, so this
 * registry needs representatives only for post-seed classes learned by
 * workers, and its global ID sequence begins after that dense seed.
 */
class sharedCanonicalIdRegistry
{
public:
    static constexpr std::size_t shardCount = 64;

    explicit sharedCanonicalIdRegistry(std::size_t firstUnusedId):
        nextId(firstUnusedId)
    {
        if (
            firstUnusedId >
            static_cast<std::size_t>(std::numeric_limits<int>::max())
        ) throw std::length_error("canonical ID space is exhausted");
    }

    sharedCanonicalIdRegistry(const sharedCanonicalIdRegistry &) = delete;
    sharedCanonicalIdRegistry &operator=(
        const sharedCanonicalIdRegistry &
    ) = delete;

    template<typename EqualRepresentative>
    int findOrInsert(
        std::size_t invariantHash,
        std::vector<std::uint64_t> words,
        EqualRepresentative &&equalRepresentative
    )
    {
        shard &selected = shards[shardIndex(invariantHash)];
        std::size_t checked = 0;
        while (true)
        {
            const representative *candidate = nullptr;
            {
                std::lock_guard lock(selected.mutex);
                auto &bucket = selected.buckets[invariantHash];
                if (checked == bucket.size())
                {
                    bucket.push_back({
                        std::move(words),
                        unknownCanonicalId
                    });
                    std::uint64_t allocated = nextId.load(
                        std::memory_order_relaxed
                    );
                    while (true)
                    {
                        if (
                            allocated > static_cast<std::uint64_t>(
                                std::numeric_limits<int>::max()
                            )
                        )
                        {
                            bucket.pop_back();
                            throw std::length_error(
                                "canonical ID space is exhausted"
                            );
                        }
                        if (nextId.compare_exchange_weak(
                            allocated,
                            allocated + 1,
                            std::memory_order_relaxed
                        )) break;
                    }
                    bucket.back().canonicalId = static_cast<int>(allocated);
                    return static_cast<int>(allocated);
                }
                candidate = std::addressof(bucket[checked++]);
            }

            // deque element addresses remain stable as concurrent insertions
            // append new immutable representatives to this bucket.
            if (equalRepresentative(candidate->words))
                return candidate->canonicalId;
        }
    }

private:
    struct representative
    {
        std::vector<std::uint64_t> words;
        int canonicalId;
    };

    struct alignas(64) shard
    {
        std::mutex mutex;
        std::unordered_map<
            std::size_t,
            std::deque<representative>
        > buckets;
    };

    std::array<shard, shardCount> shards;
    std::atomic<std::uint64_t> nextId;

    static std::size_t shardIndex(std::size_t hash) noexcept
    {
        static_assert((shardCount & (shardCount - 1)) == 0);
        return hash & (shardCount - 1);
    }
};

// Canonical maps, interners, and miss-path scratch are worker-local in OpenMP
// builds and calculation-local otherwise.
inline ASSEMBLYCPP_SEARCH_LOCAL canonicalisationGraphWorkspace
    canonicalisationGraphScratch;

inline ASSEMBLYCPP_SEARCH_LOCAL sharedCanonicalIdRegistry
    *sharedCanonicalRegistry = nullptr;

void prepareCanonicalisationGraph(
    const molGraph &source,
    const std::vector<MoleculeEdge> &edgeList
)
{
    canonicalisationGraphScratch.configure(source, edgeList);
}

#ifdef ASSEMBLYCPP_LIBRARY_BUILD
ASSEMBLYCPP_SEARCH_LOCAL std::unordered_map<
    graphHash,
    IntegerPair,
    graphHashHasher
> graphHashMap;
#else
ASSEMBLYCPP_SEARCH_LOCAL std::unordered_map<graphHash, IntegerPair>
    graphHashMap;
#endif

// The producer's canonical classes are immutable after DAG construction.
// Workers consult this shared base and insert only post-seed classes into the
// thread-local graphHashMap delta above.
inline ASSEMBLYCPP_SEARCH_LOCAL const decltype(graphHashMap)
    *sharedGraphHashSeed = nullptr;

/** Complete lazy key state once, before any worker can observe the seed. */
void freezeGraphHashSeed(const decltype(graphHashMap) &seed)
{
    for (const auto &entry : seed) entry.first.prepareForSharing();
}

void bindGraphHashSeed(const decltype(graphHashMap) &seed)
{
    if (!graphHashMap.empty())
    {
        throw std::logic_error(
            "graph hash seed requires an empty worker delta"
        );
    }
    sharedGraphHashSeed = &seed;
}

void clearGraphHashDelta() noexcept
{
    graphHashMap.clear();
    sharedGraphHashSeed = nullptr;
}

[[nodiscard]] std::size_t graphHashClassCount()
{
    const std::size_t seedSize = sharedGraphHashSeed == nullptr
        ? 0 : sharedGraphHashSeed->size();
    if (
        seedSize >
            std::numeric_limits<std::size_t>::max() - graphHashMap.size()
    ) throw std::length_error("canonical class count exceeds capacity");
    return seedSize + graphHashMap.size();
}

std::vector<std::uint64_t> serializeCanonicalMask(const EdgeMask &mask)
{
    std::vector<std::uint64_t> words(EdgeMask::activeWordCount());
    for (std::size_t word = 0; word < words.size(); ++word)
        words[word] = mask.activeWord(word);
    return words;
}

/** Keep the allocation-heavy miss path out of the cache-hit instruction body. */
ASSEMBLYCPP_NOINLINE int canoniseCacheMiss(EdgeMask &mask)
{
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.canonicalisationMaskCacheMisses;
#endif

    bool isCyclic = false;
    const molGraph &molecule = searchTargetMolecule();
    const vector<MoleculeEdge> &edgeList = searchUniverseEdgeList();
    const flatCanonGraph &graph = canonicalisationGraphScratch.build(
        molecule,
        edgeList,
        mask,
        isCyclic
    );
    graphHash candidate(graph, isCyclic);

    const IntegerPair *canonicalEntry = nullptr;
    if (sharedGraphHashSeed != nullptr)
    {
        bool canMatchSeed = true;
        if (
            sharedTreeCanonInterner != nullptr &&
            !candidate.treeHash.empty()
        )
        {
            const treeCanonNodeId seededNodeCount =
                static_cast<treeCanonNodeId>(
                    sharedTreeCanonInterner->size()
                );
            canMatchSeed =
                candidate.treeHash.first <= seededNodeCount &&
                candidate.treeHash.second <= seededNodeCount;
        }
        if (canMatchSeed)
        {
            const auto seeded = sharedGraphHashSeed->find(candidate);
            if (seeded != sharedGraphHashSeed->end())
                canonicalEntry = std::addressof(seeded->second);
        }
    }
    auto localEntry = graphHashMap.end();
    bool inserted = false;
    if (canonicalEntry == nullptr)
    {
        auto localInsertion = graphHashMap.try_emplace(
            std::move(candidate),
            IntegerPair{unknownCanonicalId, 0}
        );
        localEntry = localInsertion.first;
        inserted = localInsertion.second;
        if (inserted)
        {
            try
            {
                int canonicalId;
                if (sharedCanonicalRegistry == nullptr)
                {
                    const std::size_t classCount = graphHashClassCount();
                    const std::size_t canonicalIndex = classCount - 1;
                    if (
                        canonicalIndex > static_cast<std::size_t>(
                            std::numeric_limits<int>::max()
                        )
                    )
                    {
                        throw std::length_error(
                            "canonical ID space is exhausted"
                        );
                    }
                    canonicalId = static_cast<int>(canonicalIndex);
                }
                else
                {
                    const std::size_t invariantHash =
                        sharedCanonicalInvariantHash(graph);
                    const graphHash &localCandidate = localEntry->first;
                    canonicalId = sharedCanonicalRegistry->findOrInsert(
                        invariantHash,
                        serializeCanonicalMask(mask),
                        [&localCandidate, &molecule, &edgeList](
                            const std::vector<std::uint64_t>
                                &representativeWords
                        )
                        {
                            EdgeMask representative = EdgeMask::fromActiveWords(
                                representativeWords.data()
                            );
                            bool representativeIsCyclic = false;
                            const flatCanonGraph &representativeGraph =
                                canonicalisationGraphScratch.build(
                                    molecule,
                                    edgeList,
                                    representative,
                                    representativeIsCyclic
                                );
                            graphHash representativeHash(
                                representativeGraph,
                                representativeIsCyclic
                            );
                            return localCandidate == representativeHash;
                        }
                    );
                }
                localEntry->second = IntegerPair{canonicalId, 1};
            }
            catch (...)
            {
                graphHashMap.erase(localEntry);
                throw;
            }
        }
        canonicalEntry = std::addressof(localEntry->second);
    }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
    {
        if (inserted) ++searchTelemetry.counters.canonicalClassInsertions;
        else ++searchTelemetry.counters.canonicalClassReuses;
    }
#endif
    // Seed hit counts were used only while constructing the DAG. Keep the
    // published seed immutable; post-seed reuse counts remain worker-local.
    if (!inserted && localEntry != graphHashMap.end())
        ++localEntry->second.second;

    const IntegerPair canonical = *canonicalEntry;
    bitsetHashTable.emplace(mask, canonical);
    return canonical.first;
}

/**
 * @brief Returns unique hash val for subgraph. See Seet et al. section 4.3 Enumeration
 *
 * @param mask Boolean edgelist to be canonised
 * @return int canonical value
 */
int canonise(EdgeMask &mask)
{
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.canonicalisationCalls;
#endif
    const auto cached = bitsetHashTable.find(mask);
    if (cached != bitsetHashTable.end())
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.canonicalisationMaskCacheHits;
#endif
        return cached->second.first;
    }
    return canoniseCacheMiss(mask);
}
