enum class initialDagInsertionResult
{
    existing,
    retained,
    limitReached
};

struct initialDagInsertion
{
    initialDagInsertionResult result;
    int index = -1;
};

/**
 * @brief Insert one initial-DAG state if it is new and fits the budget.
 *
 * The level map is both the state store and the uniqueness index. A rejected
 * insertion is erased by iterator so its mask is not hashed a second time.
 */
initialDagInsertion tryRetainInitialDagMask(
    initialDagLevel &level,
    const EdgeMask &mask,
    size_t &retainedStateCount
)
{
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.retainedMaskAttempts;
#endif
    const size_t proposedIndex = level.size();
    auto [insertion, inserted] = level.try_emplace(mask);
    if (!inserted)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.duplicateMaskAttempts;
#endif
        return {initialDagInsertionResult::existing};
    }

    const size_t limit = ENUM_MAX > 0 ? static_cast<size_t>(ENUM_MAX) : 0;
    if (retainedStateCount >= limit)
    {
        level.erase(insertion);
        enumerationLimitReached = true;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.rejectedMasks;
#endif
        return {initialDagInsertionResult::limitReached};
    }
    if (proposedIndex > static_cast<size_t>(numeric_limits<int>::max()))
    {
        level.erase(insertion);
        throw length_error("initial DAG level exceeds index capacity");
    }
    ++retainedStateCount;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.retainedMasks;
#endif
    insertion->second.index = static_cast<int>(proposedIndex);
    return {initialDagInsertionResult::retained, insertion->second.index};
}

/**
 * @brief Struct for storing a potential duplicate
 */
struct potentialDuplicate
{
    /// @brief mask representing edge list of potential duplicate
    EdgeMask mask;
    /// @brief the index from the canonise function and index of the fragment
    int idx, fragment;
    potentialDuplicate() = default;

    potentialDuplicate(EdgeMask _mask, int _fragment, int _idx):
        mask(std::move(_mask)), idx(_idx), fragment(_fragment) {}
};

/**
 * @brief Struct for storing a potential duplicate during the initial enumeration before the construction of the DAG
 *
 */
struct initialPotentialDuplicate : potentialDuplicate
{
    /// mask representing presence of specific atoms in the potential duplicate
    AtomMask atomMask = 0;
    /// mask representing the edge list of the parent fragment
    EdgeMask fragMask = 0;

    /**
     * @brief Construct a new potential Duplicate object
     *
     * @param x edge to be set
     * @param _fragMask Boolean edgelist of the fragment the duplicate is part of
     * @param _fragment Index of the fragment in its assembly state
     */
    initialPotentialDuplicate(int x, EdgeMask &_fragMask, size_t _fragment)
    {
        fragMask = _fragMask;
        fragment = _fragment;
        mask.set(x);
        atomMask.set(univEdgeList[x].a);
        atomMask.set(univEdgeList[x].b);
    }

    initialPotentialDuplicate(
        const initialPotentialDuplicate &parent,
        size_t edge,
        size_t atomA,
        size_t atomB
    ):
        potentialDuplicate(
            parent.mask.withBitSet(edge),
            parent.fragment,
            parent.idx
        ),
        atomMask(parent.atomMask.withBitSet(atomA)),
        fragMask(parent.fragMask)
    {
        atomMask.set(atomB);
    }

    /**
     * @brief Generate potential matches originating from this fragment and update the DAG
     *
     * @param q Potential duplicates which are isomorphic to this.mask
     * @param retainedStateCount Number of unique masks currently held in tempDag
     * @param tempDag Temporary DAG populated with generated children
     */
    bool generateDAG(vector<initialPotentialDuplicate> &q, size_t &retainedStateCount,
    vector<initialDagLevel> &tempDag)
    {
        const size_t childLevelIndex = mask.count();
        initialDagLevel &childLevel = tempDag[childLevelIndex];
        vector<dagTransition> *transitions = nullptr;
        for (size_t i = 0; i < univEdgeList.size(); i++)
        {
            if (searchShouldStop()) return false;
            if (!mask[i] && fragMask[i])
            {
                const size_t atomA = univEdgeList[i].a;
                const size_t atomB = univEdgeList[i].b;
                if (atomMask[atomA] || atomMask[atomB])
                {
                    EdgeMask tempMask = mask.withBitSet(i);
                    const initialDagInsertion insertion =
                        tryRetainInitialDagMask(
                            childLevel,
                            tempMask,
                            retainedStateCount
                        );
                    // The initial DAG is a first-discovery forest. An existing
                    // state already has its one retained parent transition.
                    if (
                        insertion.result ==
                        initialDagInsertionResult::existing
                    )
                        continue;
                    if (
                        insertion.result ==
                        initialDagInsertionResult::limitReached
                    )
                        return false;

                    initialPotentialDuplicate g(*this, i, atomA, atomB);
                    q.push_back(std::move(g));
                    if (transitions == nullptr)
                    {
                        transitions = &tempDag[childLevelIndex - 1]
                            .at(mask)
                            .transitions;
                    }
                    transitions->emplace_back(insertion.index, i);
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
    const EdgeMask &first, &second;
    /// @brief frag1 and frag2 index first and second; maxFragSize is their maximum size
    int frag1, frag2, maxFragSize;
    validMatchings(
        const EdgeMask &_first,
        const EdgeMask &_second,
        int _frag1,
        int _frag2,
        int _maxFragSize
    ):
    first(_first), second(_second), frag1(_frag1), frag2(_frag2), maxFragSize(_maxFragSize){}
};

template<typename potentialDuplicate>
struct duplicateSet
{
    /// @brief bitset count of the duplicates in the duplicate set
    size_t size;
    /// @brief fragment masks from the assembly states
    vector<EdgeMask> maskList;
    /// @brief list of potential duplicates
    vector<potentialDuplicate> list;
    duplicateSet(size_t _size, size_t fragments)
    {
        size = _size; maskList.resize(fragments);
    }

    /**
     * @brief Insert a potential duplicate into the list
     */
    void insert(potentialDuplicate m)
    {
        maskList[m.fragment] |= m.mask;
        list.push_back(std::move(m));
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
     * @brief Visit pairable duplicates in reverse lexicographic order
     *
     * This is the order previously produced by generating every matching and
     * then consuming that vector backwards. Only the matching currently being
     * visited is materialised.
     *
     * @param visitor Called for each valid matching; false stops iteration
     * @return true if every matching was visited
     * @return false if the visitor stopped iteration or the search should stop
     */
    template<typename Visitor>
    bool visitMatchingsInReverse(Visitor &&visitor)
    {
        if (list.size() < 2) return true;

        for (size_t firstIndex = list.size() - 1; firstIndex > 0;)
        {
            --firstIndex;
            if (searchShouldStop()) return false;
            int frag = list[firstIndex].fragment;
            for (size_t secondIndex = list.size();
                 secondIndex > firstIndex + 1;)
            {
                --secondIndex;
                if (searchShouldStop()) return false;
                if (
                    frag == list[secondIndex].fragment &&
                    !list[firstIndex].mask.disjoint(list[secondIndex].mask)
                ) continue;

                validMatchings matching(
                    list[firstIndex].mask,
                    list[secondIndex].mask,
                    frag,
                    list[secondIndex].fragment,
                    size
                );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                if (searchTelemetryEnabled) [[unlikely]]
                    ++searchTelemetry.counters.matchingVisits;
#endif
                if (!visitor(matching)) return false;
            }
        }
        return true;
    }

    /**
     * @brief Check whether occurrences belong to more than one fragment
     */
    bool occurrencesSpanMultipleFragments()
    {
        bool foundFragment = false;
        for (const EdgeMask &mask : maskList)
        {
            if (searchShouldStop()) return false;
            if (mask == 0) continue;
            if (foundFragment) return true;
            foundFragment = true;
        }
        return false;
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
     * @param retainedStateCount Number of unique masks currently held in tempDag
     * @param tempDag the temporary DAG
     * @return true if any valid matchings exist
     * @return false otherwise
     */
    bool dagPopulator(vector<initialPotentialDuplicate> &q, 
    size_t &retainedStateCount,
    vector<initialDagLevel> &tempDag)
    {
        bool output = 0;
        const bool allAlive = occurrencesSpanMultipleFragments();
        if (searchShouldStop()) return output;

        auto populateDAG = [&](initialPotentialDuplicate &duplicate)
        {
            if (!duplicate.generateDAG(q, retainedStateCount, tempDag))
                return false;
            output = 1;
            return true;
        };
        if (allAlive)
        {
            for (initialPotentialDuplicate &duplicate : list)
            {
                if (searchShouldStop()) return output;
                if (!populateDAG(duplicate)) return output;
            }
            return output;
        }

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
                    if (list[i].mask.disjoint(list[j].mask))
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
                if (!populateDAG(list[i])) return output;
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
bool dagGenerate(potentialDuplicate &d, map<int, dagDuplicateSet> &stmap, EdgeMask &fragment,
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
                potentialDuplicate child(
                    d.mask.withBitSet(transition.addedEdge),
                    d.fragment,
                    d.idx
                );
                child.idx = transition.childIndex;
                auto entry = stmap.try_emplace(
                    dn.ix,
                    size + 1,
                    frags
                ).first;
                entry->second.insert(std::move(child));
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
    vector<EdgeMask> &takenMasks, vector<EdgeMask> &stateMasks, int ordinal, bool &overweight, bool last)
    {
        bool output = 0;
        const bool allAlive = ds.occurrencesSpanMultipleFragments();
        if (searchShouldStop()) return output;

        auto generateFromDuplicate = [&](potentialDuplicate &duplicate)
        {
            const int frag = duplicate.fragment;
            takenMasks[frag] |= duplicate.mask;
            ds.dead = 0;
            if (!last)
            {
                overweight |= dagGenerate(
                    duplicate,
                    stmap,
                    stateMasks[frag],
                    ds.size,
                    ordinal,
                    ds.maskList.size()
                );
                if (searchShouldStop()) return false;
                output = 1;
            }
            return true;
        };
        if (allAlive)
        {
            for (potentialDuplicate &duplicate : ds.list)
            {
                if (searchShouldStop()) return output;
                if (!generateFromDuplicate(duplicate)) return output;
            }
            return output;
        }

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
                    if (ds.list[i].mask.disjoint(ds.list[j].mask))
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
                if (!generateFromDuplicate(ds.list[i])) return output;
            }
        }
        return output;
    }
