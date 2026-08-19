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
template<>
struct std::hash<graphHash>
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
                std::hash<unsigned char>{}(gh.treeHash.centralBond)
            );
            return result;
        }
    }
};

std::unordered_map<graphHash, pii> graphHashMap;

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
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.canonicalisationMaskCacheMisses;
#endif

    bool isCyclic;
    molGraph mg = constructFromEdgeList(
        targetMolecule,
        univEdgeList,
        mask,
        isCyclic
    );
    graphHash candidate(mg, isCyclic);

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
