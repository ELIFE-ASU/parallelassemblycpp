/**
 * @brief Hashes a molecular graph
 */
struct graphHash
{
    /// graph to hash expressed as a bitset of edges
    standardBitset mask;
    /// hashes stored as floats
    vector<float> hashes;
    /// if the graph is acyclic, the tree hash function is used, and the output stored here
    string treeHash;

    graphHash() = default;

    /**
     * @brief Calculate hash for subgraph unordered_map using BFS approach
     *
     * @param mg molGraph to be hashed
     * @param _depth Depth of the BFS
     */
    void calcHash(molGraph &mg, int _depth)
    {   
        const double depthFactor = 0.33;
        int depth = min(_depth, HASH_DEPTH_MAX);
        vector<double> dhashes(hashes.size(), 0.0);
        for (size_t i = 0; i < mg.mg.size(); i++)
        {
            const string &atype1 = mg.mg[i].type;
            auto atomType = atypeHash.try_emplace(
                atype1,
                static_cast<int>((atypeHash.size() + 1) * 5)
            ).first;
            dhashes[i] = atomType->second / depthFactor;
            for (size_t j = 0; j < mg.degree(i); j++)
            {
                short btype = mg.btype(i, j);
                const string &atype = mg.mg[mg.elem(i, j)].type;
                auto neighbourType = atypeHash.try_emplace(
                    atype,
                    static_cast<int>((atypeHash.size() + 1) * 5)
                ).first;
                dhashes[i] += neighbourType->second + btype;
            }
        }
        for (int k = 1; k < depth; k++)
        {
            vector<double> oldHashes = dhashes;
            for (size_t i = 0; i < mg.mg.size(); i++)
            {
                for (size_t j = 0; j < mg.degree(i); j++)
                {
                    short btype = mg.btype(i, j);
                    dhashes[i] += (oldHashes[mg.elem(i, j)] + btype)*k*depthFactor;
                }
            }
        }
        for (size_t i = 0; i < hashes.size(); i++)
        {
            hashes[i] = static_cast<float>(dhashes[i]);
        }
    }
    
    /**
     * @brief Construct a new graph Hash object
     *
     * @param mg molGraph to be hashed
     * @param depth Depth of the BFS
     * @param isCyclic Is the molecule cyclic
     * @param _mask Boolean edgelist of the molGraph
     */
    graphHash(molGraph &mg, int depth, bool isCyclic, standardBitset &_mask)
    {
        mask = _mask;
        if (isCyclic)
        {
            hashes.resize(mg.mg.size(), 0);
            calcHash(mg, depth);
        }
        else treeHash = centroidTreeCanon(mg, 0);
    }

    /**
     * @brief Check isomorphism between two graphs. Uses tree isomorphism if acyclic, else uses vf2 subgraph isomorphism
     *
     * @param g2 other graph to be compared
     * @return true if graphs are isomorphic
     * @return false otherwise
     */
    bool operator==(const graphHash &g2) const
    {
        if (treeHash.length() != g2.treeHash.length()) return false;
        if (treeHash.length() == 0)
        {
            molGraphBoost g1mg = edgelistToBoost(targetMolecule, univEdgeList, this->mask), 
            g2mg = edgelistToBoost(targetMolecule, univEdgeList, g2.mask);
            return vf2GraphIso(g1mg, g2mg);
        }
        else
        {
            return treeHash == g2.treeHash;
        }
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
        if (gh.treeHash.length() == 0)
        {
            vector<float> sortedHashes = gh.hashes;
            sort(sortedHashes.begin(), sortedHashes.end());
            size_t res = 17;
            for (const float value : sortedHashes)
            {
                int k = int(value * 1024);
                res = res * 31 + hash<int>()(k);
            }
            return res;
        }
        else
        {
            hash<string> hasher;
            return hasher(gh.treeHash);
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
int canonise(standardBitset &mask)
{
    const auto cached = bitsetHashTable.find(mask);
    if (cached != bitsetHashTable.end()) return cached->second.first;

    bool isCyclic;
    molGraph mg = constructFromEdgeList(
        targetMolecule,
        univEdgeList,
        mask,
        isCyclic
    );
    graphHash candidate(mg, mg.mg.size(), isCyclic, mask);

    const int nextCanonicalIndex = static_cast<int>(graphHashMap.size());
    auto [graphEntry, inserted] = graphHashMap.try_emplace(
        std::move(candidate),
        pii{nextCanonicalIndex, 1}
    );
    if (!inserted) graphEntry->second.second++;

    const pii canonicalEntry = graphEntry->second;
    bitsetHashTable.emplace(mask, canonicalEntry);
    return canonicalEntry.first;
}
