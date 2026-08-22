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

    explicit atom(string _type): type(std::move(_type)){}

};

/** Exact key for edges equivalent during unique-edge preprocessing. */
struct bondClassKey
{
    string firstAtomType;
    string secondAtomType;
    short bondType;

    bool operator==(const bondClassKey &) const = default;
};

struct bondClassKeyHash
{
    static void combine(size_t &seed, size_t value)
    {
        seed ^= value + static_cast<size_t>(0x9e3779b97f4a7c15ULL)
            + (seed << 6) + (seed >> 2);
    }

    size_t operator()(const bondClassKey &key) const
    {
        size_t result = std::hash<string>{}(key.firstAtomType);
        combine(result, std::hash<short>{}(key.bondType));
        combine(result, std::hash<string>{}(key.secondAtomType));
        return result;
    }
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
    void addAtom(string _type)
    {
        mg.emplace_back(std::move(_type));
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
    size_t degree(size_t x) const
    {
        return mg[x].list.size();
    }

    short elem(size_t a, size_t b) const
    {
        return mg[a].list[b].n;
    }

    /**
     * @brief Get atom type for index i
     */
    const string &atype(size_t i) const {return mg[i].type;}

    /**
     * @brief Get bond type as short
     * @param a Index of first atom/node
     * @param b Index of bond
     */
    short btypeS(size_t a, size_t b) const
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
                        btypeS(source, bondIndex)
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
    vector<edgeL> writeEdgeList() const
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
    void writeEdgeList(
        std::unordered_map<
            bondClassKey,
            pair<int, edgeL>,
            bondClassKeyHash
        > &ht
    ) const
    {
        for (size_t i = 0; i < mg.size(); i++)
        {
            for (size_t j = 0; j < degree(i); j++)
            {
                short k = elem(i, j);
                if (i < static_cast<size_t>(k))
                {
                    const string &sourceType = atype(i);
                    const string &targetType = atype(k);
                    bondClassKey key = sourceType < targetType
                        ? bondClassKey{sourceType, targetType, btypeS(i, j)}
                        : bondClassKey{targetType, sourceType, btypeS(i, j)};
                    short source = static_cast<short>(i);
                    short bondIndex = static_cast<short>(j);
                    edgeL t(source, k, bondIndex);
                    auto [entry, inserted] = ht.try_emplace(
                        std::move(key),
                        1,
                        t
                    );
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
    /**
     * @brief For compensating for disjoint fragments in the JAI
     * 
     * @return int number of disjoint fragments
     */
    int disjointFragments() const
    {
        vb visited(mg.size(), 0);
        vector<size_t> pending;
        pending.reserve(mg.size());
        int count = 0;
        for (size_t i = 0; i < mg.size(); i++)
        {
            if (visited[i]) continue;
            count++;
            visited[i] = true;
            pending.push_back(i);
            while (!pending.empty())
            {
                const size_t vertex = pending.back();
                pending.pop_back();
                for (const bond &edge : mg[vertex].list)
                {
                    const size_t neighbour = static_cast<size_t>(edge.n);
                    if (visited[neighbour]) continue;
                    visited[neighbour] = true;
                    pending.push_back(neighbour);
                }
            }
        }
        return count;
    }
};

/**
 * @brief Preprocesses the graph by removing all unique edges for the pathway algorithm
 * @param mg The input molGraph
 * @param writeback Edges removed during preprocessing
 * @return molGraph (the final output)
 */
molGraph preprocessWriteback(const molGraph &mg, vector<edgeL> &writeback)
{
    std::unordered_map<
        bondClassKey,
        pair<int, edgeL>,
        bondClassKeyHash
    > ht;
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
    std::sort(
        v.begin(),
        v.end(),
        [](const edgeL &left, const edgeL &right)
        {
            if (left.a != right.a) return left.a < right.a;
            if (left.b != right.b) return left.b < right.b;
            return left.c < right.c;
        }
    );
    out.negativeEdgeCollapse(v);
    writeback = v;
    return out;
}

/// Global variable for the molGraph before and after preprocessing
molGraph originalMolecule, targetMolecule;

struct cachedResidualDecomposition
{
    bool isIdentity = false;
    int identityCanonicalId = unknownCanonicalId;
    vector<assemblyFragment> components;

    void appendTo(
        uint64_t maskWord,
        int activeEdgeCount,
        vector<assemblyFragment> &output
    ) const
    {
        if (isIdentity)
        {
            output.emplace_back(
                EdgeMask(maskWord),
                activeEdgeCount,
                identityCanonicalId,
                true
            );
        }
        else
        {
            output.insert(output.end(), components.begin(), components.end());
        }
    }

    void appendTo(
        const EdgeMask &mask,
        int activeEdgeCount,
        vector<assemblyFragment> &output
    ) const
    {
        if (isIdentity)
        {
            output.emplace_back(
                mask,
                activeEdgeCount,
                identityCanonicalId,
                true
            );
        }
        else
        {
            output.insert(output.end(), components.begin(), components.end());
        }
    }
};

struct lowResidualDecompositionCacheEntry
{
    uint64_t key = 0;
    uint64_t generation = 0;
    bool occupied = false;
    cachedResidualDecomposition decomposition;
};

struct wideResidualDecompositionCacheEntry
{
    uint64_t generation = 0;
    cachedResidualDecomposition decomposition;
};

enum class residualCanonicalIdCacheKind : unsigned char
{
    low,
    wide
};

/** @brief Deferred cache-ID writeback for one raw child decomposition. */
struct residualCanonicalIdBinding
{
    residualCanonicalIdCacheKind kind;
    size_t entryIndex;
    uint64_t generation;
    size_t outputStart;
};

struct ufdsMaskWorkspace
{
    // Cache molecules large enough to amortize lookup. Low-word masks retain
    // their specialised scalar table; wider masks use compact active-word keys.
    static constexpr size_t decompositionCacheEntryLimit = 4096;
    static constexpr size_t decompositionCacheComponentLimit = 16384;
    static constexpr size_t decompositionCacheMinimumMoleculeEdges = 31;
    static constexpr size_t decompositionCacheMaximumMoleculeEdges =
        EdgeMask::capacity();
    static constexpr size_t decompositionSeenBitCount = 1 << 18;
    static constexpr size_t decompositionSeenWordCount =
        decompositionSeenBitCount / numeric_limits<uint64_t>::digits;
    static constexpr size_t wideDecompositionSeenWordCount =
        decompositionCacheEntryLimit / numeric_limits<uint64_t>::digits;
    static constexpr size_t wideDecompositionUnprovenProbeLimit = 6144;
    static_assert(std::has_single_bit(decompositionCacheEntryLimit));
    static_assert(std::has_single_bit(decompositionSeenBitCount));
    static_assert(
        decompositionCacheEntryLimit < numeric_limits<uint16_t>::max()
    );

    ufdsSplit sets;
    vector<EdgeMask> components;
    vi boundTotals;
    vector<lowResidualDecompositionCacheEntry> lowDecompositionCache;
    vector<wideResidualDecompositionCacheEntry> wideDecompositionCache;
    vector<uint16_t> wideDecompositionCacheSlots;
    vector<uint64_t> wideDecompositionCacheKeys;
    vector<uint64_t> decompositionSeenBits;
    unique_ptr<uint64_t[]> wideDecompositionSeenFingerprints;
    vector<uint64_t> wideDecompositionSeenOccupied;
    vector<residualCanonicalIdBinding> residualCanonicalIdBindings;
    size_t edgeCount;
    size_t decompositionCacheKeyWordCount;
    bool reuseResidualDecompositions;
    size_t lowDecompositionCacheEntries = 0;
    size_t wideDecompositionCacheEntries = 0;
    size_t decompositionCacheComponentUnits = 0;
    uint64_t decompositionCacheGeneration = 0;
    bool decompositionCacheDisabled = false;
    size_t wideDecompositionUnprovenProbes = 0;
    bool wideDecompositionCacheProvenUseful = false;
#ifdef FRAGMENT_CACHE_STATS
    size_t decompositionCacheHits = 0;
    size_t decompositionCacheMisses = 0;
    size_t decompositionCacheAdmissions = 0;
    size_t decompositionCacheFirstOccurrences = 0;
    size_t decompositionCacheIdentityAdmissions = 0;
    size_t decompositionCacheEmptyAdmissions = 0;
    size_t wideDecompositionFingerprintProbes = 0;
    size_t wideDecompositionFirstRepeatProbe = 0;
    size_t wideDecompositionFirstHitProbe = 0;
#endif
    // Non-empty only when the processed molecule is one uniformly labelled
    // path. Values map universe edge indices to positions along that path.
    vector<int> homogeneousPathEdgePositions;

    ufdsMaskWorkspace(size_t atomCount, size_t _edgeCount):
        edgeCount(_edgeCount),
        decompositionCacheKeyWordCount(
            (_edgeCount + numeric_limits<uint64_t>::digits - 1) /
            numeric_limits<uint64_t>::digits
        ),
        reuseResidualDecompositions(decompositionCacheEligible(_edgeCount))
    {
        sets.elements.resize(atomCount);
        sets.extraVals.reserve(_edgeCount);
        components.reserve(atomCount);
        boundTotals.reserve(_edgeCount);
        residualCanonicalIdBindings.reserve(4);
        if (reuseResidualDecompositions)
        {
            if (!usesWideDecompositionCache())
                decompositionSeenBits.resize(decompositionSeenWordCount);
        }
    }

    static bool decompositionCacheEligible(size_t moleculeEdgeCount)
    {
        return
            moleculeEdgeCount >= decompositionCacheMinimumMoleculeEdges &&
            moleculeEdgeCount <= decompositionCacheMaximumMoleculeEdges;
    }

    bool usesWideDecompositionCache() const
    {
        return edgeCount > numeric_limits<uint64_t>::digits;
    }

    uint64_t nextDecompositionCacheGeneration()
    {
        ++decompositionCacheGeneration;
        if (decompositionCacheGeneration == 0) ++decompositionCacheGeneration;
        return decompositionCacheGeneration;
    }

    [[gnu::always_inline]] void beginFragmentation()
    {
        if (!residualCanonicalIdBindings.empty()) [[unlikely]]
            residualCanonicalIdBindings.clear();
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
            << " entries="
            << lowDecompositionCacheEntries + wideDecompositionCacheEntries
            << " units=" << decompositionCacheComponentUnits
            << " wide-probes=" << wideDecompositionFingerprintProbes
            << " wide-first-repeat=" << wideDecompositionFirstRepeatProbe
            << " wide-first-hit=" << wideDecompositionFirstHitProbe
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

    uint64_t wideDecompositionFingerprint(const EdgeMask &mask) const
    {
        constexpr uint64_t hashConstant = 0x9e3779b97f4a7c15ULL;
        uint64_t fingerprint =
            hashConstant ^ static_cast<uint64_t>(decompositionCacheKeyWordCount);
        for (size_t i = 0; i < decompositionCacheKeyWordCount; i++)
        {
            fingerprint ^= mask.activeWord(i) + hashConstant +
                (fingerprint << 6) + (fingerprint >> 2);
        }
        fingerprint ^= fingerprint >> 32;
        fingerprint *= 0xd6e8feb86659fd93ULL;
        fingerprint ^= fingerprint >> 32;
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

    bool wideDecompositionSeenBefore(uint64_t fingerprint)
    {
        if (wideDecompositionSeenFingerprints == nullptr)
        {
            wideDecompositionSeenFingerprints =
                make_unique_for_overwrite<uint64_t[]>(
                    decompositionCacheEntryLimit
                );
            wideDecompositionSeenOccupied.resize(
                wideDecompositionSeenWordCount
            );
        }
        constexpr size_t entryMask = decompositionCacheEntryLimit - 1;
        constexpr size_t wordShift = 6;
        const size_t index = fingerprint & entryMask;
        const uint64_t occupiedMask = uint64_t{1} << (index & 63);
        uint64_t &occupiedWord =
            wideDecompositionSeenOccupied[index >> wordShift];
        const bool seen =
            (occupiedWord & occupiedMask) != 0 &&
            wideDecompositionSeenFingerprints[index] == fingerprint;
        wideDecompositionSeenFingerprints[index] = fingerprint;
        occupiedWord |= occupiedMask;
        return seen;
    }

    bool wideDecompositionKeyEquals(
        size_t entryIndex,
        const EdgeMask &mask
    ) const
    {
        const size_t keyStart =
            entryIndex * decompositionCacheKeyWordCount;
        for (size_t i = 0; i < decompositionCacheKeyWordCount; i++)
        {
            if (wideDecompositionCacheKeys[keyStart + i] != mask.activeWord(i))
            {
                return false;
            }
        }
        return true;
    }

    void assignWideDecompositionKey(
        size_t entryIndex,
        const EdgeMask &mask
    )
    {
        const size_t keyStart =
            entryIndex * decompositionCacheKeyWordCount;
        for (size_t i = 0; i < decompositionCacheKeyWordCount; i++)
        {
            wideDecompositionCacheKeys[keyStart + i] = mask.activeWord(i);
        }
    }

    static bool needsCanonicalIdWriteback(
        const cachedResidualDecomposition &decomposition
    )
    {
        if (decomposition.isIdentity)
            return decomposition.identityCanonicalId < 0;
        // A decomposition is admitted before any of its fresh components are
        // canonicalised, then all IDs are written back together.
        return
            !decomposition.components.empty() &&
            decomposition.components.front().canonicalId < 0;
    }

    void bindLowCanonicalIds(size_t entryIndex, size_t outputStart)
    {
        const lowResidualDecompositionCacheEntry &entry =
            lowDecompositionCache[entryIndex];
        if (!needsCanonicalIdWriteback(entry.decomposition)) return;
        residualCanonicalIdBindings.push_back({
            residualCanonicalIdCacheKind::low,
            entryIndex,
            entry.generation,
            outputStart
        });
    }

    void bindWideCanonicalIds(size_t entryIndex, size_t outputStart)
    {
        const wideResidualDecompositionCacheEntry &entry =
            wideDecompositionCache[entryIndex];
        if (!needsCanonicalIdWriteback(entry.decomposition)) return;
        residualCanonicalIdBindings.push_back({
            residualCanonicalIdCacheKind::wide,
            entryIndex,
            entry.generation,
            outputStart
        });
    }

    /** Cache IDs resolved only after the raw child survives its bounds. */
    [[gnu::always_inline]] void cacheCanonicalIds(
        const vector<assemblyFragment> &output
    )
    {
        if (residualCanonicalIdBindings.empty()) [[likely]] return;
        cacheCanonicalIdsSlow(output);
    }

    [[gnu::noinline]] void cacheCanonicalIdsSlow(
        const vector<assemblyFragment> &output
    )
    {
        for (const residualCanonicalIdBinding &binding :
             residualCanonicalIdBindings)
        {
            cachedResidualDecomposition *decomposition = nullptr;
            if (binding.kind == residualCanonicalIdCacheKind::low)
            {
                if (binding.entryIndex >= lowDecompositionCache.size())
                    continue;
                lowResidualDecompositionCacheEntry &entry =
                    lowDecompositionCache[binding.entryIndex];
                if (
                    !entry.occupied ||
                    entry.generation != binding.generation
                ) continue;
                decomposition = &entry.decomposition;
            }
            else
            {
                if (binding.entryIndex >= wideDecompositionCache.size())
                    continue;
                wideResidualDecompositionCacheEntry &entry =
                    wideDecompositionCache[binding.entryIndex];
                if (entry.generation != binding.generation) continue;
                decomposition = &entry.decomposition;
            }

            const size_t componentCount = decomposition->isIdentity
                ? 1
                : decomposition->components.size();
            if (
                binding.outputStart > output.size() ||
                componentCount > output.size() - binding.outputStart
            )
            {
                throw logic_error("residual canonical-ID binding is invalid");
            }
            if (decomposition->isIdentity)
            {
                decomposition->identityCanonicalId =
                    output[binding.outputStart].canonicalId;
                continue;
            }
            for (size_t i = 0; i < componentCount; i++)
            {
                decomposition->components[i].canonicalId =
                    output[binding.outputStart + i].canonicalId;
            }
        }
        residualCanonicalIdBindings.clear();
    }

    void disableDecompositionCache()
    {
        residualCanonicalIdBindings.clear();
        vector<lowResidualDecompositionCacheEntry>().swap(
            lowDecompositionCache
        );
        vector<wideResidualDecompositionCacheEntry>().swap(
            wideDecompositionCache
        );
        vector<uint16_t>().swap(wideDecompositionCacheSlots);
        vector<uint64_t>().swap(wideDecompositionCacheKeys);
        wideDecompositionSeenFingerprints.reset();
        vector<uint64_t>().swap(wideDecompositionSeenOccupied);
        lowDecompositionCacheEntries = 0;
        wideDecompositionCacheEntries = 0;
        decompositionCacheComponentUnits = 0;
        decompositionCacheDisabled = true;
    }

    bool assignCachedComponents(
        vector<assemblyFragment> &stored,
        const vector<assemblyFragment> &output,
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
        vector<assemblyFragment> replacement(
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
 * @param fragmentList Connected residual fragments returned
 * @param workspace Reusable disjoint-set and component buffers. Its component
 * buffer must not alias fragmentList.
 */
void ufdsMaskConstructWithoutCacheWithWorkspace(
    const EdgeMask &mask,
    vector<assemblyFragment> &fragmentList,
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
    u.splitWithBuffers(fragmentList, workspace.components);
}

/**
 * @brief Append a residual mask's connected components, reusing prior results
 *
 * The cache belongs to the per-calculation workspace because an edge mask is
 * meaningful only for the molecule that produced univEdgeList.
 */
void ufdsMaskConstructWithLowCacheWithWorkspace(
    const EdgeMask &mask,
    vector<assemblyFragment> &fragmentList,
    ufdsMaskWorkspace &workspace
)
{
    const uint64_t lowMaskWord =
        bitsetLowWordBelow(mask, workspace.edgeCount);
    const size_t activeEdgeCount = std::popcount(lowMaskWord);
    auto constructWithoutCache = [&]() {
        ufdsMaskConstructWithoutCacheWithWorkspace(
            mask,
            fragmentList,
            workspace,
            lowMaskWord
        );
    };
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

    const size_t index = fingerprint &
        (ufdsMaskWorkspace::decompositionCacheEntryLimit - 1);
    const cachedResidualDecomposition *cachedDecomposition = nullptr;
    if (!workspace.decompositionCacheDisabled)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheLookups;
#endif
        if (!workspace.lowDecompositionCache.empty())
        {
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
            const size_t outputStart = fragmentList.size();
            cachedDecomposition->appendTo(
                lowMaskWord,
                static_cast<int>(activeEdgeCount),
                fragmentList
            );
            workspace.bindLowCanonicalIds(index, outputStart);
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

    const size_t outputStart = fragmentList.size();
    constructWithoutCache();

    if (
        workspace.decompositionCacheDisabled
    )
    {
        return;
    }

    const size_t componentCount = fragmentList.size() - outputStart;
    const bool isIdentity =
        componentCount == 1 &&
        fragmentList[outputStart].mask == EdgeMask(lowMaskWord);
    const size_t storedComponentCount = isIdentity ? 0 : componentCount;

    try
    {
        if (workspace.lowDecompositionCache.empty())
        {
            workspace.lowDecompositionCache.resize(
                ufdsMaskWorkspace::decompositionCacheEntryLimit
            );
        }
        lowResidualDecompositionCacheEntry &entry =
            workspace.lowDecompositionCache[index];
        vector<assemblyFragment> &stored = entry.decomposition.components;
        if (!workspace.assignCachedComponents(
            stored,
            fragmentList,
            outputStart,
            storedComponentCount
        ))
        {
            return;
        }
        entry.decomposition.isIdentity = isIdentity;
        entry.decomposition.identityCanonicalId = isIdentity
            ? fragmentList[outputStart].canonicalId
            : unknownCanonicalId;
        entry.key = lowMaskWord;
        entry.generation = workspace.nextDecompositionCacheGeneration();
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
        workspace.disableDecompositionCache();
    }
}

/**
 * @brief Append a wide residual decomposition using exact active-word keys
 *
 * A bounded fingerprint-tag table filters one-off residuals without the
 * saturation behaviour of the low-word Bloom filter. Fingerprints select a
 * direct-mapped slot, but a hit is accepted only after every active key word
 * matches.
 */
void ufdsMaskConstructWithWideCacheWithWorkspace(
    const EdgeMask &mask,
    vector<assemblyFragment> &fragmentList,
    ufdsMaskWorkspace &workspace
)
{
    const size_t activeEdgeCount = mask.count();
    auto constructWithoutCache = [&]() {
        ufdsMaskConstructWithoutCacheWithWorkspace(
            mask,
            fragmentList,
            workspace,
            0
        );
    };
    if (activeEdgeCount < 2)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheSmallResidualBypasses;
#endif
        return;
    }

    if (
        !workspace.decompositionCacheDisabled &&
        !workspace.wideDecompositionCacheProvenUseful &&
        workspace.wideDecompositionUnprovenProbes >=
            ufdsMaskWorkspace::wideDecompositionUnprovenProbeLimit
    )
    {
        workspace.disableDecompositionCache();
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
        workspace.wideDecompositionFingerprint(mask);
    if (!workspace.wideDecompositionCacheProvenUseful)
    {
        workspace.wideDecompositionUnprovenProbes++;
    }
    bool fingerprintSeen = false;
    try
    {
        fingerprintSeen =
            workspace.wideDecompositionSeenBefore(fingerprint);
    }
    catch (const bad_alloc &)
    {
        workspace.disableDecompositionCache();
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheRuntimeDisabledBypasses;
#endif
        constructWithoutCache();
        return;
    }
#ifdef FRAGMENT_CACHE_STATS
    workspace.wideDecompositionFingerprintProbes++;
    if (
        fingerprintSeen &&
        workspace.wideDecompositionFirstRepeatProbe == 0
    )
    {
        workspace.wideDecompositionFirstRepeatProbe =
            workspace.wideDecompositionFingerprintProbes;
    }
#endif
    if (!fingerprintSeen)
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
    size_t matchedEntryIndex = numeric_limits<size_t>::max();
    const size_t index = fingerprint &
        (ufdsMaskWorkspace::decompositionCacheEntryLimit - 1);
    if (!workspace.decompositionCacheDisabled)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.residualCacheLookups;
#endif
        if (!workspace.wideDecompositionCacheSlots.empty())
        {
            const uint16_t entryIndex =
                workspace.wideDecompositionCacheSlots[index];
            if (
                entryIndex != numeric_limits<uint16_t>::max() &&
                workspace.wideDecompositionKeyEquals(entryIndex, mask)
            )
            {
                matchedEntryIndex = entryIndex;
                cachedDecomposition =
                    &workspace.wideDecompositionCache[entryIndex].decomposition;
            }
        }
        if (cachedDecomposition != nullptr)
        {
            workspace.wideDecompositionCacheProvenUseful = true;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            if (searchTelemetryEnabled) [[unlikely]]
                ++searchTelemetry.counters.residualCacheHits;
#endif
#ifdef FRAGMENT_CACHE_STATS
            workspace.decompositionCacheHits++;
            if (workspace.wideDecompositionFirstHitProbe == 0)
            {
                workspace.wideDecompositionFirstHitProbe =
                    workspace.wideDecompositionFingerprintProbes;
            }
#endif
            const size_t outputStart = fragmentList.size();
            cachedDecomposition->appendTo(
                mask,
                static_cast<int>(activeEdgeCount),
                fragmentList
            );
            workspace.bindWideCanonicalIds(matchedEntryIndex, outputStart);
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

    const size_t outputStart = fragmentList.size();
    constructWithoutCache();

    if (workspace.decompositionCacheDisabled) return;

    const size_t componentCount = fragmentList.size() - outputStart;
    const bool isIdentity =
        componentCount == 1 && fragmentList[outputStart].mask == mask;
    const size_t storedComponentCount = isIdentity ? 0 : componentCount;

    try
    {
        if (workspace.wideDecompositionCacheSlots.empty())
        {
            workspace.wideDecompositionCacheSlots.resize(
                ufdsMaskWorkspace::decompositionCacheEntryLimit,
                numeric_limits<uint16_t>::max()
            );
        }
        uint16_t &cacheSlot = workspace.wideDecompositionCacheSlots[index];
        size_t entryIndex = cacheSlot;
        const bool newEntry = cacheSlot == numeric_limits<uint16_t>::max();
        if (newEntry)
        {
            entryIndex = workspace.wideDecompositionCache.size();
            workspace.wideDecompositionCache.emplace_back();
            workspace.wideDecompositionCacheKeys.resize(
                (entryIndex + 1) * workspace.decompositionCacheKeyWordCount
            );
        }
        wideResidualDecompositionCacheEntry &entry =
            workspace.wideDecompositionCache[entryIndex];
        vector<assemblyFragment> &stored = entry.decomposition.components;
        if (!workspace.assignCachedComponents(
            stored,
            fragmentList,
            outputStart,
            storedComponentCount
        ))
        {
            if (newEntry)
            {
                workspace.wideDecompositionCache.pop_back();
                workspace.wideDecompositionCacheKeys.resize(
                    entryIndex * workspace.decompositionCacheKeyWordCount
                );
            }
            return;
        }
        entry.decomposition.isIdentity = isIdentity;
        entry.decomposition.identityCanonicalId = isIdentity
            ? fragmentList[outputStart].canonicalId
            : unknownCanonicalId;
        entry.generation = workspace.nextDecompositionCacheGeneration();
        workspace.assignWideDecompositionKey(entryIndex, mask);
        if (newEntry)
        {
            cacheSlot = static_cast<uint16_t>(entryIndex);
            workspace.wideDecompositionCacheEntries++;
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
        workspace.disableDecompositionCache();
    }
}

/**
 * @brief Append a residual decomposition, using reuse only where it pays
 */
inline void ufdsMaskConstructWithWorkspace(
    const EdgeMask &mask,
    vector<assemblyFragment> &fragmentList,
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
        if (workspace.decompositionCacheDisabled)
        {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            if (searchTelemetryEnabled) [[unlikely]]
            {
                if (mask.count() < 2)
                {
                    ++searchTelemetry.counters
                        .residualCacheSmallResidualBypasses;
                }
                else
                {
                    ++searchTelemetry.counters
                        .residualCacheRuntimeDisabledBypasses;
                }
            }
#endif
            ufdsMaskConstructWithoutCacheWithWorkspace(
                mask,
                fragmentList,
                workspace,
                workspace.usesWideDecompositionCache()
                    ? 0
                    : bitsetLowWordBelow(mask, workspace.edgeCount)
            );
            return;
        }
        if (workspace.usesWideDecompositionCache())
        {
            ufdsMaskConstructWithWideCacheWithWorkspace(
                mask,
                fragmentList,
                workspace
            );
        }
        else
        {
            ufdsMaskConstructWithLowCacheWithWorkspace(
                mask,
                fragmentList,
                workspace
            );
        }
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
        fragmentList,
        workspace,
        workspace.edgeCount <= numeric_limits<uint64_t>::digits
            ? bitsetLowWordBelow(mask, workspace.edgeCount)
            : 0
    );
}
