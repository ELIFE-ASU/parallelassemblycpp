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
     * @param mg molGraph to be hashed
     * @param isCyclic Is the molecule cyclic
     */
    graphHash(molGraph &mg, bool isCyclic)
    {
        if (isCyclic) cyclicHash = canonicaliseCyclicGraph(mg);
        else
        {
            treeHash = centroidTreeCanon(mg, 0);
            if (treeHash.empty()) cyclicHash = canonicaliseWholeGraph(mg);
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
    const edgeL *configuredEdgeData = nullptr;

    /** Must be called after changing source atoms or the universe edge list. */
    void configure(const molGraph &source, const std::vector<edgeL> &edgeList)
    {
        if (source.mg.size() > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("canonical graph has too many vertices");

        sourceLabels.clear();
        sourceLabels.reserve(source.mg.size());
        sourceLegacyX.clear();
        sourceLegacyX.reserve(source.mg.size());
        for (const atom &vertex : source.mg)
        {
            sourceLabels.push_back(internTreeCanonAtom(vertex.type));
            sourceLegacyX.push_back(vertex.type == "X");
        }

        sourceEdges.clear();
        sourceEdges.reserve(edgeList.size());
        for (const edgeL &edge : edgeList)
        {
            if (
                edge.a < 0 || edge.b < 0 || edge.c < 0 || edge.a == edge.b ||
                static_cast<std::size_t>(edge.a) >= source.mg.size() ||
                static_cast<std::size_t>(edge.b) >= source.mg.size() ||
                static_cast<std::size_t>(edge.c) >=
                    source.mg[static_cast<std::size_t>(edge.a)].list.size()
            )
            {
                throw std::logic_error("canonical source edge is invalid");
            }
            const bond &sourceBond = source.mg[edge.a].list[edge.c];
            sourceEdges.push_back({
                static_cast<std::uint32_t>(edge.a),
                static_cast<std::uint32_t>(edge.b),
                canonGraphBondType(sourceBond)
            });
        }

        localVertex.resize(source.mg.size());
        vertexEpoch.resize(source.mg.size(), 0);
        unionNodes.reserve(source.mg.size());
        degree.reserve(source.mg.size());
        selectedEdges.reserve(edgeList.size());
        graph.labels.reserve(source.mg.size());
        graph.adjacencyOffsets.reserve(source.mg.size() + 1);
        graph.adjacency.reserve(edgeList.size() * 2);
        configuredInternerGeneration = treeCanonInternerGeneration;
        configuredAtomData = source.mg.data();
        configuredEdgeData = edgeList.data();
    }

    void ensureConfigured(
        const molGraph &source,
        const std::vector<edgeL> &edgeList
    )
    {
        // This catches interner resets and storage replacement. Call configure
        // explicitly after any in-place atom or edge content mutation, even
        // when vector sizes and backing addresses are unchanged.
        if (
            configuredInternerGeneration != treeCanonInternerGeneration ||
            sourceLabels.size() != source.mg.size() ||
            sourceEdges.size() != edgeList.size() ||
            configuredAtomData != source.mg.data() ||
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
    [[gnu::noinline]] void addSelectedEdge(
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
        const std::vector<edgeL> &edgeList,
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

// Enumeration, canonical maps, and interners are already process-global and
// single-threaded; the flat miss-path storage follows the same lifecycle.
inline canonicalisationGraphWorkspace canonicalisationGraphScratch;

void prepareCanonicalisationGraph(
    const molGraph &source,
    const std::vector<edgeL> &edgeList
)
{
    canonicalisationGraphScratch.configure(source, edgeList);
}

#ifdef ASSEMBLYCPP_LIBRARY_BUILD
std::unordered_map<graphHash, pii, graphHashHasher> graphHashMap;
#else
std::unordered_map<graphHash, pii> graphHashMap;
#endif

/** Keep the allocation-heavy miss path out of the cache-hit instruction body. */
[[gnu::noinline]] int canoniseCacheMiss(EdgeMask &mask)
{
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.canonicalisationMaskCacheMisses;
#endif

    bool isCyclic = false;
    const flatCanonGraph &graph = canonicalisationGraphScratch.build(
        targetMolecule,
        univEdgeList,
        mask,
        isCyclic
    );
    graphHash candidate(graph, isCyclic);

    const int nextCanonicalIndex = static_cast<int>(graphHashMap.size());
    auto [graphEntry, inserted] = graphHashMap.try_emplace(
        std::move(candidate),
        pii{nextCanonicalIndex, 1}
    );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
    {
        if (inserted) ++searchTelemetry.counters.canonicalClassInsertions;
        else ++searchTelemetry.counters.canonicalClassReuses;
    }
#endif
    if (!inserted) graphEntry->second.second++;

    const pii canonicalEntry = graphEntry->second;
    bitsetHashTable.emplace(mask, canonicalEntry);
    return canonicalEntry.first;
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
