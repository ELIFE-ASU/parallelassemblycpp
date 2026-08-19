/**
 * @brief Bond struct for molGraph
 */
struct bond
{
    short n;
    short type;
    bond(short _n, short _type): n(_n), type(_type){}
};

/**
 * @brief Atom struct for molGraph
 */
struct atom
{
    string type;
    vector<bond> list;

    atom(string _type): type(_type){}

};

/**
 * @brief Primary graph data structure used in assemblyCpp
 */
struct molGraph
{
    /**
     * @brief Vector of atoms representing nodes in the graph.
     */
    vector<atom> mg;
    /**
     * @brief Total number of bonds (edges) in the graph.
     */
    int totalBonds = 0;

    /**
     * @brief Use this function to add atoms/nodes
     * @param _type Type is atom type/node labelling.
     */
    void addAtom(string &_type)
    {
        atom a(_type);
        mg.push_back(a);
    }

    /**
     * @brief Use this function to add bonds/edges.
     * @param a Index of first atom/node
     * @param b Index of second atom/node
     * @param type Type is bond order/edge labelling.
     */
    void addBond(int a, int b, short type)
    {
        bond b1(b, type), b2(a, type);
        mg[a].list.push_back(b1);
        mg[b].list.push_back(b2);
        totalBonds++;
    }
    
    /**
     * @brief Print the graph information to cout
     */
    void printToCout()
    {
        cout << "There are " << mg.size() << " atoms in the molecule-graph\n";
        for (size_t i = 0; i < mg.size(); i++)
        {
            cout << "Atom " << i + 1 << " is of type " << mg[i].type << " and adjacent to atoms ";
            for (size_t j = 0; j < degree(i); j++)
            {
                cout << elem(i, j) + 1 << " with bond order " << btypeS(i, j) << ", ";
            }
            cout << '\n';
        }
    }

    /**
     * @brief Get the degree (number of bonds) of atom at index x
     */
    size_t degree(int x)
    {
        return mg[x].list.size();
    }

    short elem(size_t a, size_t b)
    {
        return mg[a].list[b].n;
    }

    /**
     * @brief Get atom type for index i
     */
    string atype(size_t i) {return mg[i].type;}

    /**
     * @brief Get bond type as char
     * @param a Index of first atom/node
     * @param b Index of bond
     */
    char btype(size_t a, size_t b)
    {
        return static_cast<char>(mg[a].list[b].type);
    }

    /**
     * @brief Get bond type as short
     * @param a Index of first atom/node
     * @param b Index of bond
     */
    short btypeS(size_t a, size_t b)
    {
        return mg[a].list[b].type;
    }

private:
    /**
     * @brief Rebuild the graph while dropping removed atoms and zero-order bonds.
     */
    void rebuild(bool removeMarkedAtoms)
    {
        const size_t removed = mg.size();
        vector<size_t> reverseMap(mg.size(), removed);
        molGraph output;

        for (size_t i = 0; i < mg.size(); i++)
        {
            if (removeMarkedAtoms && mg[i].type == "COLLAPSE") continue;
            reverseMap[i] = output.mg.size();
            output.addAtom(mg[i].type);
        }

        for (size_t source = 0; source < mg.size(); source++)
        {
            if (reverseMap[source] == removed) continue;

            for (size_t bondIndex = 0; bondIndex < degree(source); bondIndex++)
            {
                const size_t target = static_cast<size_t>(elem(source, bondIndex));
                if (
                    reverseMap[target] != removed &&
                    btypeS(source, bondIndex) != 0 &&
                    reverseMap[source] < reverseMap[target]
                )
                {
                    output.addBond(
                        static_cast<int>(reverseMap[source]),
                        static_cast<int>(reverseMap[target]),
                        btype(source, bondIndex)
                    );
                }
            }
        }

        *this = std::move(output);
    }

public:

    /**
     * @brief Remove all bonds of order 0. Used in preprocessing.
     */
    void collapse()
    {
        rebuild(false);
    }

    /**
     * @brief For explicit hydrogen removal
     * 
     */
    void removeAtom(size_t i)
    {
        if (i >= mg.size()) return;
        mg[i].type = "COLLAPSE";
    }

    /**
     * @brief For explicit hydrogen removal
     * 
     */
    void removeAndCollapse()
    {
        rebuild(true);
    }

    /**
     * @brief Turns molGraph (adjacency list) into equivalent edgelist
     * @return std::vector<edgeL>
     */
    vector<edgeL> writeEdgeList()
    {
        vector<edgeL> out;
        for (size_t i = 0; i < mg.size(); i++)
        {
            for (size_t j = 0; j < degree(i); j++)
            {
                short k = elem(i, j);
                if (i < static_cast<size_t>(k))
                {
                    short source = static_cast<short>(i);
                    short bondIndex = static_cast<short>(j);
                    edgeL t(source, k, bondIndex);
                    out.push_back(t);
                }
            }
        }
        return out;
    }
    
    /**
     * @brief For preprocessing, writes edgeList as hash map to detect duplicated bonds
     */
    void writeEdgeList(std::unordered_map<string, pair<int, edgeL> > &ht)
    {
        for (size_t i = 0; i < mg.size(); i++)
        {
            for (size_t j = 0; j < degree(i); j++)
            {
                short k = elem(i, j);
                if (i < static_cast<size_t>(k))
                {
                    string is = atype(i), ks = atype(k), out;
                    if (is < ks) out = is + btype(i, j) + ks;
                    else out = ks + btype(i, j) + is;
                    short source = static_cast<short>(i);
                    short bondIndex = static_cast<short>(j);
                    edgeL t(source, k, bondIndex);
                    auto [entry, inserted] = ht.try_emplace(out, 1, t);
                    if (!inserted) entry->second.first++;
                }
            }
        }
    }

    /**
     * @brief Used in preprocessing, removes edges in edgelist
     */
    void negativeEdgeCollapse(vector<edgeL> &edgeList)
    {
        for (size_t i = 0; i < edgeList.size(); i++)
        {
            edgeL &el = edgeList[i];
            mg[el.a].list[el.c].type = 0;
        }
        collapse();
    }

    /**
     * @brief For compensating for disjoint fragments in the JAI
     * 
     */
    void disjointFragmentsR(vb & visited, int n)
    {
        if (visited[n]) return;
        visited[n] = 1;
        for (size_t i = 0; i < mg[n].list.size(); i++) disjointFragmentsR(visited, elem(n, i));
    }

    /**
     * @brief For compensating for disjoint fragments in the JAI
     * 
     * @return int number of disjoint fragments
     */
    int disjointFragments()
    {
        vb visited(mg.size(), 0); int count = 0;
        for (size_t i = 0; i < mg.size(); i++)
        {
            if (!visited[i])
            {
                count++;
                disjointFragmentsR(visited, i);
            }
        }
        return count;
    }
};

/**
 * @brief construct new molGraph from input molGraph and boolean edgelist
 *
 * @param mg Input molgraph
 * @param edgeList Corresponding edge list
 * @param mask Input boolean edgelist
 * @param isCyclic Is the graph cyclic, needed for hashing
 * @return molGraph
 */
molGraph constructFromEdgeList(molGraph &mg, vector<edgeL> &edgeList, 
    EdgeMask &mask, bool &isCyclic)
    {
        disjointSet u(mg.mg.size());
        molGraph output;
        std::unordered_map<int, int> ht;
        isCyclic = 0;
        for (size_t i = 0; i < edgeList.size(); i++)
        {
            if (mask[i] != 0)
            {
                int a = edgeList[i].a, b = edgeList[i].b, c = 0;
                if (ht.count(a) == 0 && ht.count(b) == 0)
                {
                    size_t x = ht.size();
                    ht[a] = x;
                    output.addAtom(mg.mg[a].type);
                    size_t y = ht.size();
                    ht[b] = y;
                    output.addAtom(mg.mg[b].type);
                    u.insert(x, x);
                    u.insert(y, x);
                }
                else
                {
                    if (ht.count(a) == 0)
                    {
                        size_t x = ht.size();
                        ht[a] = x;
                        output.addAtom(mg.mg[a].type);
                        u.insert(x, ht[b]);
                    }
                    else c++;
                    if (ht.count(b) == 0)
                    {
                        size_t x = ht.size();
                        ht[b] = x;
                        output.addAtom(mg.mg[b].type);
                        u.insert(x, ht[a]);
                    }
                    else c++;
                }
                int a2 = ht[a], b2 = ht[b];
                output.addBond(a2, b2, mg.btype(a, edgeList[i].c));
                if (c == 2)
                {
                    isCyclic |= u.merge(a2, b2);
                }
            }
        }
        return output;
    }

/**
 * @brief Preprocesses the graph by removing all unique edges for the pathway algorithm
 * @param mg The input molGraph
 * @param writeback Edges removed during preprocessing
 * @return molGraph (the final output)
 */
molGraph preprocessWriteback(molGraph &mg, vector<edgeL> &writeback)
{
    std::unordered_map<string, pair<int, edgeL> > ht;
    molGraph out = mg;
    mg.writeEdgeList(ht);
    vector<edgeL> v;
    for (auto it = ht.begin(); it != ht.end(); ++it)
    {
        if (it->second.first == 1)
        {
            v.push_back(it->second.second);
        }
    }
    out.negativeEdgeCollapse(v);
    writeback = v;
    return out;
}

/// Global variable for the molGraph before and after preprocessing
molGraph originalMolecule, targetMolecule;

struct cachedResidualDecomposition
{
    bool isIdentity = false;
    vector<EdgeMask> components;
    vi edgeCounts;

    void appendTo(
        uint64_t maskWord,
        int activeEdgeCount,
        vector<EdgeMask> &output,
        vi &outputEdgeCounts
    ) const
    {
        if (isIdentity)
        {
            output.emplace_back(maskWord);
            outputEdgeCounts.push_back(activeEdgeCount);
        }
        else
        {
            output.insert(output.end(), components.begin(), components.end());
            outputEdgeCounts.insert(
                outputEdgeCounts.end(),
                edgeCounts.begin(),
                edgeCounts.end()
            );
        }
    }
};

struct lowResidualDecompositionCacheEntry
{
    uint64_t key = 0;
    bool occupied = false;
    cachedResidualDecomposition decomposition;
};

struct ufdsMaskWorkspace
{
    // Cache only low-word molecules large enough to amortize the lookup.
    // Smaller and wider masks retain the direct residual-proportional path.
    static constexpr size_t decompositionCacheEntryLimit = 4096;
    static constexpr size_t decompositionCacheComponentLimit = 16384;
    static constexpr size_t decompositionCacheMinimumMoleculeEdges = 31;
    static constexpr size_t decompositionCacheMaximumMoleculeEdges =
        numeric_limits<uint64_t>::digits;
    static constexpr size_t decompositionSeenBitCount = 1 << 18;
    static constexpr size_t decompositionSeenWordCount =
        decompositionSeenBitCount / numeric_limits<uint64_t>::digits;
    static_assert(std::has_single_bit(decompositionCacheEntryLimit));
    static_assert(std::has_single_bit(decompositionSeenBitCount));

    ufdsSplit sets;
    vector<EdgeMask> components;
    vi boundTotals;
    vector<lowResidualDecompositionCacheEntry> lowDecompositionCache;
    vector<uint64_t> decompositionSeenBits;
    size_t edgeCount;
    bool reuseResidualDecompositions;
    size_t lowDecompositionCacheEntries = 0;
    size_t decompositionCacheComponentUnits = 0;
    bool decompositionCacheDisabled = false;
#ifdef FRAGMENT_CACHE_STATS
    size_t decompositionCacheHits = 0;
    size_t decompositionCacheMisses = 0;
    size_t decompositionCacheAdmissions = 0;
    size_t decompositionCacheFirstOccurrences = 0;
    size_t decompositionCacheIdentityAdmissions = 0;
    size_t decompositionCacheEmptyAdmissions = 0;
#endif

    ufdsMaskWorkspace(size_t atomCount, size_t _edgeCount):
        edgeCount(_edgeCount),
        reuseResidualDecompositions(decompositionCacheEligible(_edgeCount))
    {
        sets.elements.resize(atomCount);
        sets.extraVals.reserve(_edgeCount);
        components.reserve(atomCount);
        boundTotals.reserve(_edgeCount);
        if (reuseResidualDecompositions)
        {
            decompositionSeenBits.resize(decompositionSeenWordCount);
        }
    }

    static bool decompositionCacheEligible(size_t moleculeEdgeCount)
    {
        return
            moleculeEdgeCount >= decompositionCacheMinimumMoleculeEdges &&
            moleculeEdgeCount <= decompositionCacheMaximumMoleculeEdges;
    }

#ifdef FRAGMENT_CACHE_STATS
    ~ufdsMaskWorkspace()
    {
        cerr
            << "fragment-cache hits=" << decompositionCacheHits
            << " misses=" << decompositionCacheMisses
            << " admissions=" << decompositionCacheAdmissions
            << " first=" << decompositionCacheFirstOccurrences
            << " identity=" << decompositionCacheIdentityAdmissions
            << " empty=" << decompositionCacheEmptyAdmissions
            << " entries=" << lowDecompositionCacheEntries
            << " units=" << decompositionCacheComponentUnits
            << '\n';
    }
#endif

    static uint64_t mixDecompositionFingerprint(uint64_t fingerprint)
    {
        fingerprint ^= fingerprint >> 30;
        fingerprint *= 0xbf58476d1ce4e5b9ULL;
        fingerprint ^= fingerprint >> 27;
        fingerprint *= 0x94d049bb133111ebULL;
        fingerprint ^= fingerprint >> 31;
        return fingerprint;
    }

    bool decompositionSeenBefore(uint64_t fingerprint)
    {
        constexpr size_t bitMask = decompositionSeenBitCount - 1;
        constexpr size_t wordShift = 6;
        const size_t firstBit = fingerprint & bitMask;
        const size_t secondBit = (fingerprint >> 32) & bitMask;
        const uint64_t firstMask = uint64_t{1} << (firstBit & 63);
        const uint64_t secondMask = uint64_t{1} << (secondBit & 63);
        uint64_t &firstWord = decompositionSeenBits[firstBit >> wordShift];
        uint64_t &secondWord = decompositionSeenBits[secondBit >> wordShift];
        const bool seen =
            (firstWord & firstMask) != 0 &&
            (secondWord & secondMask) != 0;
        firstWord |= firstMask;
        secondWord |= secondMask;
        return seen;
    }

    bool assignCachedComponents(
        vector<EdgeMask> &stored,
        const vector<EdgeMask> &output,
        size_t outputStart,
        size_t storedComponentCount
    )
    {
        const size_t oldCapacity = stored.capacity();
        if (storedComponentCount == 0)
        {
            stored.clear();
            return true;
        }
        if (storedComponentCount <= oldCapacity)
        {
            stored.resize(storedComponentCount);
            copy(
                output.begin() + outputStart,
                output.end(),
                stored.begin()
            );
            return true;
        }

        const size_t unitsWithoutStored =
            decompositionCacheComponentUnits - oldCapacity;
        if (
            storedComponentCount >
            decompositionCacheComponentLimit - unitsWithoutStored
        )
        {
            return false;
        }
        vector<EdgeMask> replacement(
            output.begin() + outputStart,
            output.end()
        );
        const size_t newUnits = unitsWithoutStored + replacement.capacity();
        if (newUnits > decompositionCacheComponentLimit) return false;
        stored.swap(replacement);
        decompositionCacheComponentUnits = newUnits;
        return true;
    }
};

/**
 * @brief Calls the disjoint-set data structure for the fragmentation function. See Seet et al. section 4.5 for details
 *
 * @param mask Target bitset as input
 * @param maskList List of disjoint bitsets returned
 * @param workspace Reusable disjoint-set and component buffers. Its component
 * buffer must not alias maskList.
 */
void ufdsMaskConstructWithoutCacheWithWorkspace(
    const EdgeMask &mask,
    vector<EdgeMask> &maskList,
    vi &edgeCounts,
    ufdsMaskWorkspace &workspace,
    uint64_t lowMaskWord
)
{
    vector<edgeL> &edgeList = univEdgeList;
    ufdsSplit &u = workspace.sets;
    auto processEdge = [&](size_t i) {
        int a = edgeList[i].a, b = edgeList[i].b;
        const bool containsA = u.contains(a), containsB = u.contains(b);
        if (!containsA && !containsB)
        {
            u.doubleInsert(b, a, i);
        }
        else if (!containsA)
        {
            u.insert(a, b, i);
        }
        else if (!containsB)
        {
            u.insert(b, a, i);
        }
        else
        {
            u.merge(a, b, i);
        }
    };

    size_t firstEdge = edgeList.size();
    bool initialised = false;
    auto visitEdge = [&](size_t i) {
        if (firstEdge == edgeList.size())
        {
            firstEdge = i;
            return;
        }
        if (!initialised)
        {
            u.reset();
            processEdge(firstEdge);
            initialised = true;
        }
        processEdge(i);
    };
    if (edgeList.size() <= numeric_limits<uint64_t>::digits)
    {
        uint64_t word = lowMaskWord;
        while (word != 0)
        {
            const size_t index = std::countr_zero(word);
            visitEdge(index);
            word &= word - 1;
        }
    }
    else forEachSetBitWithWideLimit(mask, edgeList.size(), visitEdge);
    if (!initialised) return;
    u.splitWithBuffers(maskList, edgeCounts, workspace.components);
}

/**
 * @brief Append a residual mask's connected components, reusing prior results
 *
 * The cache belongs to the per-calculation workspace because an edge mask is
 * meaningful only for the molecule that produced univEdgeList.
 */
void ufdsMaskConstructWithCacheWithWorkspace(
    const EdgeMask &mask,
    vector<EdgeMask> &maskList,
    vi &edgeCounts,
    ufdsMaskWorkspace &workspace
)
{
    const uint64_t lowMaskWord =
        bitsetLowWordBelow(mask, workspace.edgeCount);
    auto constructWithoutCache = [&]() {
        ufdsMaskConstructWithoutCacheWithWorkspace(
            mask,
            maskList,
            edgeCounts,
            workspace,
            lowMaskWord
        );
    };
    const size_t activeEdgeCount = std::popcount(lowMaskWord);
    if (activeEdgeCount < 2)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheSmallResidualBypasses;
#endif
        return;
    }

    if (
        workspace.decompositionCacheDisabled ||
        activeEdgeCount < 4
    )
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
        {
            if (workspace.decompositionCacheDisabled)
            {
                ++searchTelemetry.counters
                    .residualCacheRuntimeDisabledBypasses;
            }
            else
            {
                ++searchTelemetry.counters
                    .residualCacheSmallResidualBypasses;
            }
        }
#endif
        constructWithoutCache();
        return;
    }

    const uint64_t fingerprint =
        ufdsMaskWorkspace::mixDecompositionFingerprint(
            lowMaskWord
        );
    if (!workspace.decompositionSeenBefore(fingerprint))
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheFirstOccurrenceBypasses;
#endif
#ifdef FRAGMENT_CACHE_STATS
        workspace.decompositionCacheFirstOccurrences++;
#endif
        constructWithoutCache();
        return;
    }

    const cachedResidualDecomposition *cachedDecomposition = nullptr;
    if (!workspace.decompositionCacheDisabled)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheLookups;
#endif
        if (!workspace.lowDecompositionCache.empty())
        {
            const size_t index = fingerprint &
                (ufdsMaskWorkspace::decompositionCacheEntryLimit - 1);
            const lowResidualDecompositionCacheEntry &cached =
                workspace.lowDecompositionCache[index];
            if (cached.occupied && cached.key == lowMaskWord)
            {
                cachedDecomposition = &cached.decomposition;
            }
        }
        if (cachedDecomposition != nullptr)
        {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            if (searchTelemetryEnabled) [[unlikely]]
                ++searchTelemetry.counters.residualCacheHits;
#endif
#ifdef FRAGMENT_CACHE_STATS
            workspace.decompositionCacheHits++;
#endif
            cachedDecomposition->appendTo(
                lowMaskWord,
                static_cast<int>(activeEdgeCount),
                maskList,
                edgeCounts
            );
            return;
        }
#ifdef FRAGMENT_CACHE_STATS
        workspace.decompositionCacheMisses++;
#endif
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheMisses;
#endif
    }

    const size_t outputStart = maskList.size();
    constructWithoutCache();

    if (
        workspace.decompositionCacheDisabled
    )
    {
        return;
    }

    const size_t componentCount = maskList.size() - outputStart;
    const bool isIdentity =
        componentCount == 1 &&
        maskList[outputStart] == EdgeMask(lowMaskWord);
    const size_t storedComponentCount = isIdentity ? 0 : componentCount;

    try
    {
        if (workspace.lowDecompositionCache.empty())
        {
            workspace.lowDecompositionCache.resize(
                ufdsMaskWorkspace::decompositionCacheEntryLimit
            );
        }
        const size_t index = fingerprint &
            (ufdsMaskWorkspace::decompositionCacheEntryLimit - 1);
        lowResidualDecompositionCacheEntry &entry =
            workspace.lowDecompositionCache[index];
        vector<EdgeMask> &stored = entry.decomposition.components;
        if (!workspace.assignCachedComponents(
            stored,
            maskList,
            outputStart,
            storedComponentCount
        ))
        {
            return;
        }
        entry.decomposition.edgeCounts.assign(
            edgeCounts.begin() + outputStart,
            edgeCounts.end()
        );
        entry.decomposition.isIdentity = isIdentity;
        entry.key = lowMaskWord;
        if (!entry.occupied)
        {
            entry.occupied = true;
            workspace.lowDecompositionCacheEntries++;
        }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheAdmissions;
#endif
#ifdef FRAGMENT_CACHE_STATS
        workspace.decompositionCacheAdmissions++;
        workspace.decompositionCacheIdentityAdmissions += isIdentity;
        workspace.decompositionCacheEmptyAdmissions += componentCount == 0;
#endif
    }
    catch (const bad_alloc &)
    {
        vector<lowResidualDecompositionCacheEntry>().swap(
            workspace.lowDecompositionCache
        );
        workspace.lowDecompositionCacheEntries = 0;
        workspace.decompositionCacheComponentUnits = 0;
        workspace.decompositionCacheDisabled = true;
    }
}

/**
 * @brief Append a residual decomposition, using reuse only where it pays
 */
inline void ufdsMaskConstructWithWorkspace(
    const EdgeMask &mask,
    vector<EdgeMask> &maskList,
    vi &edgeCounts,
    ufdsMaskWorkspace &workspace
)
{
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.residualDecompositionRequests;
#endif
    if (workspace.reuseResidualDecompositions)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheEligibleRequests;
#endif
        ufdsMaskConstructWithCacheWithWorkspace(
            mask,
            maskList,
            edgeCounts,
            workspace
        );
        return;
    }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
    {
        if (
            workspace.edgeCount <
            ufdsMaskWorkspace::decompositionCacheMinimumMoleculeEdges
        )
        {
            ++searchTelemetry.counters.residualCacheSmallMoleculeBypasses;
        }
        else
        {
            ++searchTelemetry.counters.residualCacheWideMoleculeBypasses;
        }
    }
#endif
    ufdsMaskConstructWithoutCacheWithWorkspace(
        mask,
        maskList,
        edgeCounts,
        workspace,
        workspace.edgeCount <= numeric_limits<uint64_t>::digits
            ? bitsetLowWordBelow(mask, workspace.edgeCount)
            : 0
    );
}
