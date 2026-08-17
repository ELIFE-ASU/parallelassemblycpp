/**
 * @brief Reserve one unique state in the initial enumeration DAG.
 *
 * Existing masks do not consume the budget twice. A new mask is added only
 * when doing so keeps the number of retained states at or below ENUM_MAX.
 */
bool reserveInitialDagMask(
    std::unordered_set<standardBitset> &maskMap,
    const standardBitset &mask
)
{
    if (maskMap.count(mask) != 0) return true;

    const size_t limit = ENUM_MAX > 0 ? static_cast<size_t>(ENUM_MAX) : 0;
    if (maskMap.size() >= limit)
    {
        enumerationLimitReached = true;
        return false;
    }
    maskMap.insert(mask);
    return true;
}

/**
 * @brief Struct for storing a potential duplicate
 */
struct potentialDuplicate
{
    /// @brief mask representing edge list of potential duplicate
    standardBitset mask;
    /// @brief the index from the canonise function and index of the fragment
    int idx, fragment;
    potentialDuplicate() = default;

    potentialDuplicate(standardBitset &_mask, int _fragment, int _idx):
        mask(_mask), idx(_idx), fragment(_fragment) {}
};

/**
 * @brief Struct for storing a potential duplicate during the initial enumeration before the construction of the DAG
 *
 */
struct initialPotentialDuplicate : potentialDuplicate
{
    /// mask representing presence of specific atoms in the potential duplicate
    standardBitset atomMask = 0;
    /// mask representing the edge list of the parent fragment
    standardBitset fragMask = 0;

    /**
     * @brief Construct a new potential Duplicate object
     *
     * @param x edge to be set
     * @param _fragMask Boolean edgelist of the fragment the duplicate is part of
     * @param _fragment Index of the fragment in its assembly state
     */
    initialPotentialDuplicate(int x, standardBitset &_fragMask, size_t _fragment)
    {
        fragMask = _fragMask;
        fragment = _fragment;
        mask.set(x);
        atomMask.set(univEdgeList[x].a);
        atomMask.set(univEdgeList[x].b);
    }

    /**
     * @brief Generate potential matches originating from this fragment and update the DAG
     *
     * @param q Potential duplicates which are isomorphic to this.mask
     * @param maskMap Hash table of boolean edgelists of all fragments taken before to avoid repetition
     * @param tempDag Temporary DAG populated with generated children
     */
    bool generateDAG(vector<initialPotentialDuplicate> &q, std::unordered_set<standardBitset> &maskMap,
    vector<std::unordered_map<standardBitset, pair<int, vector<standardBitset> > > > &tempDag)
    {
        vector<edgeL> &edgeList = univEdgeList;
        for (size_t i = 0; i < edgeList.size(); i++)
        {
            if (searchShouldStop()) return false;
            if ((mask[i] == 0) && (fragMask[i] != 0))
            {
                standardBitset temp1, temp2; temp1.set(edgeList[i].a), 
                temp2.set(edgeList[i].b);
                if ((temp1 & atomMask) == temp1 || (temp2 & atomMask) == temp2)
                {
                    standardBitset tempMask = mask; tempMask.set(i);
                    if (maskMap.count(tempMask) == 0)
                    {
                        if (!reserveInitialDagMask(maskMap, tempMask)) return false;
                        initialPotentialDuplicate g = *this;
                        g.mask.set(i);
                        g.atomMask |= (temp1 | temp2);
                        q.push_back(g);
                        vector<standardBitset> &adjList = tempDag[mask.count() - 1][mask].second;
                        pair<int, vector<standardBitset> > p;
                        tempDag[mask.count()][g.mask] = p;
                        adjList.push_back(g.mask);
                    }
                }
            }
        }
        return true;
    }
};

/**
 * @brief Struct containing pair of valid duplicates for subsequent fragmentation
 *
 */
struct validMatchings
{
    /// @brief first and second duplicate masks
    standardBitset first, second;
    /// @brief frag1 and frag2 index first and second; maxFragSize is their maximum size
    int frag1, frag2, maxFragSize;
    validMatchings(standardBitset &_first, standardBitset &_second, int _frag1, int _frag2, int _maxFragSize):
    first(_first), second(_second), frag1(_frag1), frag2(_frag2), maxFragSize(_maxFragSize){}
};

template<typename potentialDuplicate>
struct duplicateSet
{
    /// @brief bitset count of the duplicates in the duplicate set
    size_t size;
    /// @brief fragment masks from the assembly states
    vector<standardBitset> maskList;
    /// @brief list of potential duplicates
    vector<potentialDuplicate> list;
    duplicateSet(size_t _size, size_t fragments)
    {
        size = _size; maskList.resize(fragments, 0);
    }

    /**
     * @brief Insert a potential duplicate into the list
     */
    void insert(const potentialDuplicate &m)
    {
        list.push_back(m);
        maskList[m.fragment] |= m.mask;
    }

    /**
     * @brief Check if there is a valid matching in the current maskList
     *
     * @return true if there is such a matching
     * @return false otherwise
     */
    bool isValid()
    {
        int count = 0, last = 0;
        for (size_t i = 0; i < maskList.size(); i++)
        {
            if (searchShouldStop()) return false;
            if (maskList[i] != 0) {count++; last = i;}
        }
        if (count > 1) return true;
        if (maskList[last].count() < (size<<1)) return false;
        return true;
    }
    
    /**
     * @brief Generate matchings from pairable duplicates from the maskList
     * 
     * @param v the output list of valid matchings
     */
    void generateMatchings(vector<validMatchings> &v)
    {
        for (size_t i = 0; i < list.size(); i++)
        {
            if (searchShouldStop()) return;
            int frag = list[i].fragment;
            for (size_t j = i + 1; j < list.size(); j++)
            {
                if (searchShouldStop()) return;
                if (frag == list[j].fragment)
                {
                    if ((list[i].mask & list[j].mask) == 0)
                    {
                        validMatchings p(list[i].mask, list[j].mask, frag, frag, size);
                        v.push_back(p);
                    }
                }
                else
                {
                    validMatchings p(list[i].mask, list[j].mask,
                    list[i].fragment, list[j].fragment, size);
                    v.push_back(p);
                }
            }
        }
    }
};

/**
 * @brief Set of boolean edgelists which are isomorphic
 *
 */
struct initialDuplicateSet : duplicateSet<initialPotentialDuplicate>
{
    using duplicateSet::duplicateSet;
    /**
     * @brief Generate size + 1 matchings from the current set and populate the DAG during the initial enumeration
     * 
     * @param q list of potential duplicates
     * @param maskMap hash map to prevent bitsets from being added to the DAG more than once
     * @param tempDag the temporary DAG
     * @return true if any valid matchings exist
     * @return false otherwise
     */
    bool dagPopulator(vector<initialPotentialDuplicate> &q, 
    std::unordered_set<standardBitset> &maskMap, 
    vector<std::unordered_map<standardBitset, pair<int, vector<standardBitset> > > > &tempDag)
    {
        bool output = 0;
        vb alive(list.size(), 0);
        for (size_t i = 0; i < list.size(); i++)
        {
            if (searchShouldStop()) return output;
            int frag = list[i].fragment;
            for (size_t j = i + 1; j < list.size(); j++)
            {
                if (searchShouldStop()) return output;
                if (frag == list[j].fragment)
                {
                    if ((list[i].mask & list[j].mask) == 0)
                    {
                        alive[i] = 1;
                        alive[j] = 1;
                    }
                }
                else
                {
                    alive[i] = 1;
                    alive[j] = 1;
                }
            }
            if (alive[i])
            {
                if (!list[i].generateDAG(q, maskMap, tempDag)) return output;
                output = 1;
            }
        }
        return output;
    }
};

/**
 * @brief Version of initialDuplicateSet which uses the DAG to search the duplicatable subgraph space more quickly
 */
struct dagDuplicateSet : duplicateSet<potentialDuplicate>
{
    bool dead = 1;
    using duplicateSet::duplicateSet;
};

/**
 * @brief Generate the next set of duplicates from the duplicate d
 * 
 * @param d the duplicate from which the next set of duplicates is to be generated
 * @param stmap maps an integer corresponding to a unique index of each non-isomorphic graph to a set of potential duplicates
 * @param fragment the bitset corresponding to the fragment d is a part of
 * @param size the maximum allowed size of a duplicate
 * @param ordinal the maximum allowed index of a duplicate
 * @param frags the number of fragments in the assembly state
 * @return true if the canonical index of any duplicate is greater than the ordinal
 * @return false otherwise
 */
bool dagGenerate(potentialDuplicate &d, map<int, dagDuplicateSet> &stmap, standardBitset &fragment, 
    size_t size, int ordinal, size_t frags)
{
    bool overweight = 0;
    const dagNode &parent = DAG[size - 1][d.idx];
    for (const dagTransition &transition : parent.transitions)
    {
        if (searchShouldStop()) return overweight;
        if (fragment[transition.addedEdge])
        {
            dagNode &dn = DAG[size][transition.childIndex];
            if (dn.ix <= ordinal)
            {
                potentialDuplicate child = d;
                child.mask.set(transition.addedEdge);
                child.idx = transition.childIndex;
                auto entry = stmap.try_emplace(
                    dn.ix,
                    size + 1,
                    frags
                ).first;
                entry->second.insert(child);
            }
            else overweight = 1;
        }
    }
    return overweight;
}

/**
 * @brief Generate the next set of duplicates from the duplicate set ds using the function dagGenerate
 * 
 * @param ds the duplicate set from which the next set is to be generated
 * @param stmap maps an integer corresponding to a unique index of each non-isomorphic graph to a set of potential duplicates
 * @param takenMasks bitsets of all edges which could be part of a duplicate
 * @param stateMasks the bitsets of the original assembly state
 * @param ordinal the maximum allowed index of a duplicate
 * @param overweight true if the generation function has reached states which are larger than the ordinal
 * @param last true if this is to be the final iteration (previous iteration was overweight)
 * @return true if any valid duplicatable subgraphs found and not the final iteration
 * @return false 
 */
bool dagDuplicateGenerator(dagDuplicateSet &ds, map<int, dagDuplicateSet> &stmap,
    vector<standardBitset> &takenMasks, vector<standardBitset> &stateMasks, int ordinal, bool &overweight, bool last)
    {
        bool output = 0;
        vb alive(ds.list.size(), 0);
        for (size_t i = 0; i < ds.list.size(); i++)
        {
            if (searchShouldStop()) return output;
            int frag = ds.list[i].fragment;
            for (size_t j = i + 1; j < ds.list.size(); j++)
            {
                if (searchShouldStop()) return output;
                if (frag == ds.list[j].fragment)
                {
                    if ((ds.list[i].mask & ds.list[j].mask) == 0)
                    {
                        alive[i] = 1;
                        alive[j] = 1;
                    }
                }
                else
                {
                    alive[i] = 1;
                    alive[j] = 1;
                }
            }
            if (alive[i])
            {
                takenMasks[frag] |= ds.list[i].mask;
                ds.dead = 0;
                if (!last)
                {
                    overweight |= dagGenerate(ds.list[i], stmap, stateMasks[frag], ds.size,
                         ordinal, ds.maskList.size());
                    if (searchShouldStop()) return output;
                    output = 1;
                }
            }
        }
        return output;
    }
