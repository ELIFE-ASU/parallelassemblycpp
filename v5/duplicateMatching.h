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
 * The level hash is the uniqueness index while nodes are stored densely by
 * first-discovery index. A rejected insertion is erased by iterator so its
 * mask is not hashed a second time.
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
    const size_t proposedIndex = level.nodes.size();
    if (proposedIndex > static_cast<size_t>(numeric_limits<int>::max()))
        throw length_error("initial DAG level exceeds index capacity");
    auto [insertion, inserted] = level.maskIndices.try_emplace(
        mask,
        static_cast<int>(proposedIndex)
    );
    if (!inserted)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.duplicateMaskAttempts;
#endif
        return {initialDagInsertionResult::existing, insertion->second};
    }

    const size_t limit = ENUM_MAX > 0 ? static_cast<size_t>(ENUM_MAX) : 0;
    if (retainedStateCount >= limit)
    {
        level.maskIndices.erase(insertion);
        enumerationLimitReached = true;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.rejectedMasks;
#endif
        return {initialDagInsertionResult::limitReached};
    }
    level.nodes.emplace_back();
    ++retainedStateCount;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.retainedMasks;
#endif
    return {
        initialDagInsertionResult::retained,
        static_cast<int>(proposedIndex)
    };
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
    initialPotentialDuplicate(
        int x,
        EdgeMask &_fragMask,
        size_t _fragment,
        int _idx
    )
    {
        fragMask = _fragMask;
        fragment = static_cast<int>(_fragment);
        idx = _idx;
        mask.set(x);
        atomMask.set(univEdgeList[x].a);
        atomMask.set(univEdgeList[x].b);
    }

    initialPotentialDuplicate(
        const initialPotentialDuplicate &parent,
        size_t edge,
        size_t atomA,
        size_t atomB,
        int _idx
    ):
        potentialDuplicate(
            parent.mask.withBitSet(edge),
            parent.fragment,
            _idx
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
        initialDagLevel &parentLevel = tempDag[childLevelIndex - 1];
        initialDagNode *parentNode = nullptr;
        for (size_t i = 0; i < univEdgeList.size(); i++)
        {
            if (searchShouldStopPeriodically()) return false;
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

                    initialPotentialDuplicate g(
                        *this,
                        i,
                        atomA,
                        atomB,
                        insertion.index
                    );
                    q.push_back(std::move(g));
                    if (parentNode == nullptr)
                    {
                        if (
                            idx < 0 ||
                            static_cast<size_t>(idx) >= parentLevel.nodes.size()
                        )
                        {
                            throw logic_error("initial DAG parent is missing");
                        }
                        parentNode = &parentLevel.nodes[idx];
                        if (
                            parentNode->transitionOffset !=
                            unassignedDagTransitionOffset
                        )
                        {
                            throw logic_error(
                                "initial DAG parent was expanded twice"
                            );
                        }
                        if (
                            parentLevel.transitions.size() >=
                            numeric_limits<uint32_t>::max()
                        )
                        {
                            throw length_error(
                                "initial DAG transition table is too large"
                            );
                        }
                        parentNode->transitionOffset = static_cast<uint32_t>(
                            parentLevel.transitions.size()
                        );
                    }
                    if (
                        parentNode->transitionCount ==
                        numeric_limits<uint32_t>::max()
                    )
                    {
                        throw length_error(
                            "initial DAG node has too many transitions"
                        );
                    }
                    parentLevel.transitions.emplace_back(insertion.index, i);
                    ++parentNode->transitionCount;
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
    first(_first), second(_second), frag1(_frag1), frag2(_frag2),
    maxFragSize(_maxFragSize) {}
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
    template<typename Filter, typename Visitor>
    [[gnu::noinline]] bool visitMatchingsInReverse(
        Filter &&filter,
        Visitor &&visitor
    )
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
                if (searchShouldStopPeriodically()) return false;
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
                if (filter(matching, firstIndex, secondIndex)) continue;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                if (searchTelemetryEnabled) [[unlikely]]
                    ++searchTelemetry.counters.matchingVisits;
#endif
                if (!visitor(matching)) return false;
            }
        }
        return true;
    }

    template<typename Visitor>
    [[gnu::always_inline]] bool visitMatchingsInReverse(Visitor &&visitor)
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
                if (searchShouldStopPeriodically()) return false;
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
                if (searchShouldStopPeriodically()) return output;
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

/** @brief One compact canonical-class bucket in an enumeration level. */
template<typename DuplicateSetType>
struct duplicateClassEntry
{
    int canonicalId;
    DuplicateSetType duplicates;

    duplicateClassEntry(
        int _canonicalId,
        size_t duplicateSize,
        size_t fragmentCount
    ):
        canonicalId(_canonicalId),
        duplicates(duplicateSize, fragmentCount) {}
};

/**
 * @brief Compact duplicate classes, sealed in ascending canonical-ID order.
 */
template<typename DuplicateSetType>
struct duplicateClassLevel
{
    vector<duplicateClassEntry<DuplicateSetType>> classes;

    bool empty() const noexcept {return classes.empty();}
    size_t size() const noexcept {return classes.size();}

    void seal()
    {
        if (classes.size() < 2) return;
        sort(
            classes.begin(),
            classes.end(),
            [](const auto &left, const auto &right)
            {
                return left.canonicalId < right.canonicalId;
            }
        );
    }
};

using initialDuplicateClassLevel = duplicateClassLevel<initialDuplicateSet>;
using dagDuplicateClassLevel = duplicateClassLevel<dagDuplicateSet>;

/**
 * @brief Reusable dense canonical-ID lookup for one level being populated.
 *
 * The slots contain only primitive positions. Enumeration levels retain
 * ownership of their buckets while recursive child searches reuse this index.
 */
struct duplicateClassIndexWorkspace
{
    struct slot
    {
        uint32_t generation = 0;
        uint32_t position = 0;
    };

    vector<slot> slots;
    uint32_t generation = 0;

    void beginLevel()
    {
        ++generation;
        if (generation != 0) return;
        for (slot &entry : slots) entry.generation = 0;
        generation = 1;
    }

    template<typename DuplicateSetType>
    DuplicateSetType &getOrCreate(
        duplicateClassLevel<DuplicateSetType> &level,
        int canonicalId,
        size_t duplicateSize,
        size_t fragmentCount
    )
    {
        if (generation == 0) [[unlikely]]
            throw logic_error("duplicate-class level was not started");
        if (
            !level.classes.empty() &&
            level.classes.back().canonicalId == canonicalId
        ) [[likely]]
        {
            return level.classes.back().duplicates;
        }
        const size_t id = static_cast<size_t>(canonicalId);
        if (canonicalId >= 0 && id < slots.size()) [[likely]]
        {
            const slot &entry = slots[id];
            if (entry.generation == generation) [[likely]]
            {
                return level.classes[entry.position].duplicates;
            }
        }
        return create(
            level,
            canonicalId,
            duplicateSize,
            fragmentCount
        );
    }

private:
    template<typename DuplicateSetType>
    [[gnu::noinline]] DuplicateSetType &create(
        duplicateClassLevel<DuplicateSetType> &level,
        int canonicalId,
        size_t duplicateSize,
        size_t fragmentCount
    )
    {
        if (generation == 0)
            throw logic_error("duplicate-class level was not started");
        if (canonicalId < 0)
            throw logic_error("negative canonical duplicate class");
        const size_t id = static_cast<size_t>(canonicalId);
        if (id >= slots.size()) slots.resize(id + 1);
        if (level.classes.size() > numeric_limits<uint32_t>::max())
            throw length_error("duplicate-class level exceeds index capacity");

        slot &entry = slots[id];
        entry.generation = generation;
        entry.position = static_cast<uint32_t>(level.classes.size());
        level.classes.emplace_back(
            canonicalId,
            duplicateSize,
            fragmentCount
        );
        return level.classes.back().duplicates;
    }
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
bool dagGenerate(
    potentialDuplicate &d,
    dagDuplicateClassLevel &stmap,
    duplicateClassIndexWorkspace &classIndex,
    EdgeMask &fragment,
    size_t size,
    int ordinal,
    size_t frags
)
{
    bool overweight = 0;
    const dagLevel &parentLevel = DAG[size - 1];
    const dagNode &parent = parentLevel.nodes[d.idx];
    if (parent.transitionCount == 0) return overweight;
    const dagTransition *transition =
        parentLevel.transitions.data() + parent.transitionOffset;
    const dagTransition *const transitionEnd =
        transition + parent.transitionCount;
    const dagNode *const childNodes = DAG[size].nodes.data();
    for (; transition != transitionEnd; ++transition)
    {
        if (searchShouldStopPeriodically()) return overweight;
        if (fragment[transition->addedEdge])
        {
            const dagNode &dn = childNodes[transition->childIndex];
            if (dn.ix <= ordinal)
            {
                potentialDuplicate child(
                    d.mask.withBitSet(transition->addedEdge),
                    d.fragment,
                    d.idx
                );
                child.idx = transition->childIndex;
                classIndex.getOrCreate(
                    stmap,
                    dn.ix,
                    size + 1,
                    frags
                ).insert(std::move(child));
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
bool dagDuplicateGenerator(
    dagDuplicateSet &ds,
    dagDuplicateClassLevel &stmap,
    duplicateClassIndexWorkspace &classIndex,
    vector<EdgeMask> &takenMasks,
    vector<EdgeMask> &stateMasks,
    int ordinal,
    bool &overweight,
    bool last
)
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
                    classIndex,
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
                if (searchShouldStopPeriodically()) return output;
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
