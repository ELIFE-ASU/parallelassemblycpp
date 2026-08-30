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
 * @brief Compact atom-to-incident-edge index for frontier expansion.
 *
 * A full EdgeMask per atom would use quadratic auxiliary space on sparse
 * graphs. This CSR index stores the same incidence relation in O(V + E) space
 * while states carry their own eligible-edge frontier mask.
 */
struct initialIncidentEdgeIndex
{
    vector<size_t> offsets;
    vector<uint32_t> edgeIndices;
    vector<EdgeMask> singleWordMasks;

    initialIncidentEdgeIndex(): offsets(AtomMask::size() + 1, 0)
    {
        const vector<edgeL> &edgeList = searchUniverseEdgeList();
        if (
            edgeList.size() >
            static_cast<size_t>(numeric_limits<int>::max())
        )
        {
            throw length_error("edge universe exceeds frontier index capacity");
        }

        // A true per-atom edge mask is both compact and faster when the whole
        // edge domain fits in one machine word. Wider domains use the sparse
        // CSR representation below so auxiliary storage remains O(V + E).
        if (EdgeMask::activeWordCount() <= 1)
        {
            singleWordMasks.resize(AtomMask::size());
            for (size_t edgeIndex = 0;
                 edgeIndex < edgeList.size();
                 edgeIndex++)
            {
                const edgeL &edge = edgeList[edgeIndex];
                validateEndpoint(edge.a);
                validateEndpoint(edge.b);
                singleWordMasks[edge.a].set(edgeIndex);
                if (edge.b != edge.a)
                    singleWordMasks[edge.b].set(edgeIndex);
            }
            return;
        }

        for (const edgeL &edge : edgeList)
        {
            validateEndpoint(edge.a);
            validateEndpoint(edge.b);
            ++offsets[static_cast<size_t>(edge.a) + 1];
            if (edge.b != edge.a)
                ++offsets[static_cast<size_t>(edge.b) + 1];
        }
        for (size_t atom = 1; atom < offsets.size(); atom++)
            offsets[atom] += offsets[atom - 1];

        edgeIndices.resize(offsets.back());
        vector<size_t> next = offsets;
        for (size_t edgeIndex = 0; edgeIndex < edgeList.size(); edgeIndex++)
        {
            const edgeL &edge = edgeList[edgeIndex];
            const size_t atomA = static_cast<size_t>(edge.a);
            const size_t atomB = static_cast<size_t>(edge.b);
            edgeIndices[next[atomA]++] = static_cast<uint32_t>(edgeIndex);
            if (edge.b != edge.a)
                edgeIndices[next[atomB]++] = static_cast<uint32_t>(edgeIndex);
        }
    }

    [[gnu::noinline]] void addEligibleEdges(
        size_t atomA,
        size_t atomB,
        const EdgeMask &fragmentMask,
        const EdgeMask &selectedMask,
        EdgeMask &frontier
    ) const
    {
        if (atomA + 1 >= offsets.size() || atomB + 1 >= offsets.size())
            throw logic_error("frontier atom is outside the incidence index");
        if (!singleWordMasks.empty())
        {
            frontier |= singleWordMasks[atomA];
            if (atomB != atomA) frontier |= singleWordMasks[atomB];
            frontier &= fragmentMask;
            frontier &= ~selectedMask;
            return;
        }
        addEligibleEdges(atomA, fragmentMask, selectedMask, frontier);
        if (atomB != atomA)
            addEligibleEdges(atomB, fragmentMask, selectedMask, frontier);
    }

private:
    void addEligibleEdges(
        size_t atom,
        const EdgeMask &fragmentMask,
        const EdgeMask &selectedMask,
        EdgeMask &frontier
    ) const
    {
        for (size_t position = offsets[atom];
             position < offsets[atom + 1];
             position++)
        {
            const size_t edge = edgeIndices[position];
            if (fragmentMask[edge] && !selectedMask[edge]) frontier.set(edge);
        }
    }

    void validateEndpoint(short endpoint) const
    {
        if (
            endpoint < 0 ||
            static_cast<size_t>(endpoint) + 1 >= offsets.size()
        )
        {
            throw logic_error("edge endpoint is outside the atom domain");
        }
    }
};

/**
 * @brief Struct for storing a potential duplicate during the initial enumeration before the construction of the DAG
 *
 */
struct initialPotentialDuplicate : potentialDuplicate
{
    /// unselected fragment edges incident to at least one selected atom
    EdgeMask frontier = 0;

    /**
     * @brief Construct a new potential Duplicate object
     *
     * @param x edge to be set
     * @param fragmentMask Edge mask of the fragment containing the duplicate
     * @param incidentEdges Precomputed compact atom-to-edge incidence index
     * @param _fragment Index of the fragment in its assembly state
     */
    initialPotentialDuplicate(
        int x,
        const EdgeMask &fragmentMask,
        const initialIncidentEdgeIndex &incidentEdges,
        size_t _fragment,
        int _idx
    )
    {
        fragment = static_cast<int>(_fragment);
        idx = _idx;
        mask.set(x);
        const vector<edgeL> &edgeList = searchUniverseEdgeList();
        const size_t atomA = edgeList[x].a;
        const size_t atomB = edgeList[x].b;
        incidentEdges.addEligibleEdges(
            atomA,
            atomB,
            fragmentMask,
            mask,
            frontier
        );
    }

    initialPotentialDuplicate(
        const initialPotentialDuplicate &parent,
        EdgeMask childMask,
        EdgeMask childFrontier,
        int _idx
    ):
        potentialDuplicate(
            std::move(childMask),
            parent.fragment,
            _idx
        ),
        frontier(std::move(childFrontier))
    {}

    /**
     * @brief Generate potential matches originating from this fragment and update the DAG
     *
     * @param q Potential duplicates which are isomorphic to this.mask
     * @param retainedStateCount Number of unique masks currently held in tempDag
     * @param tempDag Temporary DAG populated with generated children
     */
    bool generateDAG(
        vector<initialPotentialDuplicate> &q,
        size_t &retainedStateCount,
        vector<initialDagLevel> &tempDag,
        const EdgeMask &fragmentMask,
        const initialIncidentEdgeIndex &incidentEdges
    )
    {
        const vector<edgeL> &edgeList = searchUniverseEdgeList();
        const size_t childLevelIndex = mask.count();
        initialDagLevel &childLevel = tempDag[childLevelIndex];
        initialDagLevel &parentLevel = tempDag[childLevelIndex - 1];
        initialDagNode *parentNode = nullptr;
        constexpr size_t frontierWordBits =
            numeric_limits<unsigned long long>::digits;
        for (size_t wordIndex = 0;
             wordIndex < EdgeMask::activeWordCount();
             wordIndex++)
        {
            unsigned long long frontierWord = frontier.activeWord(wordIndex);
            while (frontierWord != 0)
            {
                const size_t i = wordIndex * frontierWordBits +
                    static_cast<size_t>(std::countr_zero(frontierWord));
                frontierWord &= frontierWord - 1;
                if (searchShouldStopPeriodically()) return false;
                const size_t atomA = edgeList[i].a;
                const size_t atomB = edgeList[i].b;
                EdgeMask tempMask = mask.withBitSet(i);
                const initialDagInsertion insertion =
                    tryRetainInitialDagMask(
                        childLevel,
                        tempMask,
                        retainedStateCount
                    );
                // The initial DAG is a first-discovery forest. An existing
                // state already has its one retained parent transition.
                if (insertion.result == initialDagInsertionResult::existing)
                    continue;
                if (insertion.result == initialDagInsertionResult::limitReached)
                    return false;

                EdgeMask childFrontier = frontier;
                childFrontier.reset(i);
                incidentEdges.addEligibleEdges(
                    atomA,
                    atomB,
                    fragmentMask,
                    tempMask,
                    childFrontier
                );

                initialPotentialDuplicate g(
                    *this,
                    std::move(tempMask),
                    std::move(childFrontier),
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

template<typename PotentialDuplicate>
struct duplicateSet
{
    using occurrence_type = PotentialDuplicate;

    /// @brief bitset count of the duplicates in the duplicate set
    size_t size = 0;
    /// @brief total fragments in the assembly state (including empty CSR rows)
    size_t fragmentCount = 0;
    /// @brief stable occurrence slice owned by the enclosing class level
    span<PotentialDuplicate> list;

    /**
     * @brief Check whether the sealed class contains a valid matching
     *
     * @return true if there is such a matching
     * @return false otherwise
     */
    bool isValid()
    {
        return valid;
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
     * @brief Check whether fragment runs are dense enough for block traversal
     *
     * DAG enumeration is fragment-major and preserves occurrence order while
     * extending duplicate classes, so production lists are already grouped.
     * Unexpected interleaving, or a block layout too sparse to amortise the
     * extra dispatch, should retain legacy occurrence-pair traversal.
     */
    [[gnu::noinline]] bool hasDenseFragmentRuns() const
    {
        if (list.size() < 2) return false;

        int previousFragment = list.front().fragment;
        size_t fragmentRunCount = 1;
        for (size_t occurrence = 1; occurrence < list.size(); ++occurrence)
        {
            if (searchShouldStopPeriodically()) return false;
            const int fragment = list[occurrence].fragment;
            if (fragment < previousFragment) return false;
            fragmentRunCount += fragment != previousFragment;
            previousFragment = fragment;
        }

        // Below 16 occurrences, only a single run has enough measured reuse
        // to offset the block dispatch.
        if (list.size() < 16) return fragmentRunCount == 1;

        // For larger classes, C(n, 2) >= 4 * C(g + 1, 2)
        // simplifies to g <= (n - 2) / 2.
        return fragmentRunCount <= (list.size() - 2) / 2;
    }

    /**
     * @brief Visit Cartesian blocks formed by fragment-contiguous runs
     *
     * @pre list is grouped into nondecreasing fragment runs
     *
     * fragmentPairFilter is called exactly once per block containing a valid
     * occurrence pair and receives only the block's fragment ids and duplicate
     * size; true skips the entire block. The callback therefore cannot
     * accidentally depend on a representative occurrence's masks.
    */
    template<typename FragmentPairFilter, typename Visitor>
    [[gnu::noinline]] bool visitMatchingsByFragmentPairInReverse(
        FragmentPairFilter &&fragmentPairFilter,
        Visitor &&visitor
    )
    {
        if (list.size() < 2) return true;

        size_t firstEnd = list.size();
        while (firstEnd > 0)
        {
            if (searchShouldStop()) return false;
            const int firstFragment = list[firstEnd - 1].fragment;
            size_t firstBegin = firstEnd - 1;
            while (
                firstBegin > 0 &&
                list[firstBegin - 1].fragment == firstFragment
            )
            {
                if (searchShouldStopPeriodically()) return false;
                --firstBegin;
            }

            size_t secondEnd = list.size();
            while (secondEnd > firstBegin)
            {
                if (searchShouldStopPeriodically()) return false;
                const int secondFragment = list[secondEnd - 1].fragment;
                size_t secondBegin = secondEnd - 1;
                while (
                    secondBegin > firstBegin &&
                    list[secondBegin - 1].fragment == secondFragment
                )
                {
                    if (searchShouldStopPeriodically()) return false;
                    --secondBegin;
                }
                const bool sameFragmentBlock = secondBegin == firstBegin;
                size_t firstPositionEnd = firstEnd;
                size_t nextSecondPositionEnd = secondEnd;
                if (sameFragmentBlock)
                {
                    bool matchingExists = false;
                    for (size_t firstPosition = firstEnd - 1;
                         firstPosition > firstBegin;)
                    {
                        --firstPosition;
                        for (size_t secondPosition = firstEnd;
                             secondPosition > firstPosition + 1;)
                        {
                            --secondPosition;
                            if (searchShouldStopPeriodically()) return false;
                            if (list[firstPosition].mask.disjoint(
                                list[secondPosition].mask
                            ))
                            {
                                firstPositionEnd = firstPosition + 1;
                                nextSecondPositionEnd = secondPosition + 1;
                                matchingExists = true;
                                break;
                            }
                        }
                        if (matchingExists) break;
                    }
                    // This diagonal is the last block for the current run.
                    if (!matchingExists) break;
                }
                if (!fragmentPairFilter(
                    firstFragment,
                    secondFragment,
                    size
                ))
                {
                    for (size_t firstPosition = firstPositionEnd;
                         firstPosition > firstBegin;)
                    {
                        --firstPosition;
                        const size_t secondPositionBegin =
                            sameFragmentBlock
                                ? firstPosition + 1
                                : secondBegin;
                        for (size_t secondPosition = nextSecondPositionEnd;
                             secondPosition > secondPositionBegin;)
                        {
                            --secondPosition;
                            if (searchShouldStopPeriodically()) return false;
                            if (
                                sameFragmentBlock &&
                                !list[firstPosition].mask.disjoint(
                                    list[secondPosition].mask
                                )
                            ) continue;
                            validMatchings matching(
                                list[firstPosition].mask,
                                list[secondPosition].mask,
                                firstFragment,
                                secondFragment,
                                size
                            );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                            if (searchTelemetryEnabled) [[unlikely]]
                                ++searchTelemetry.counters.matchingVisits;
#endif
                            if (!visitor(matching)) return false;
                        }
                        nextSecondPositionEnd = secondEnd;
                    }
                }
                if (secondBegin == firstBegin) break;
                secondEnd = secondBegin;
            }
            firstEnd = firstBegin;
        }
        return true;
    }

    /**
     * @brief Check whether occurrences belong to more than one fragment
     */
    bool occurrencesSpanMultipleFragments()
    {
        return spansMultipleFragments;
    }

    void bind(
        size_t duplicateSize,
        size_t fragments,
        span<PotentialDuplicate> occurrences,
        bool isValidClass,
        bool spansFragments
    ) noexcept
    {
        size = duplicateSize;
        fragmentCount = fragments;
        list = occurrences;
        valid = isValidClass;
        spansMultipleFragments = spansFragments;
    }

private:
    bool valid = false;
    bool spansMultipleFragments = false;
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
    vector<initialDagLevel> &tempDag,
    const vector<assemblyFragment> &fragments,
    const initialIncidentEdgeIndex &incidentEdges,
    vector<uint8_t> &aliveScratch)
    {
        bool output = 0;
        const bool allAlive = occurrencesSpanMultipleFragments();
        if (searchShouldStop()) return output;

        auto populateDAG = [&](initialPotentialDuplicate &duplicate)
        {
            const bool completed = duplicate.generateDAG(
                q,
                retainedStateCount,
                tempDag,
                fragments[duplicate.fragment].mask,
                incidentEdges
            );
            // Each retained initial occurrence is expanded at most once. Its
            // frontier is not needed by matching, so release the wide mask for
            // reuse while keeping the occurrence mask and fragment identity.
            duplicate.frontier.reset();
            if (!completed) return false;
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

        aliveScratch.assign(list.size(), 0);
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
                        aliveScratch[i] = 1;
                        aliveScratch[j] = 1;
                    }
                }
                else
                {
                    aliveScratch[i] = 1;
                    aliveScratch[j] = 1;
                }
            }
            if (aliveScratch[i])
            {
                if (!populateDAG(list[i])) return output;
            }
            else list[i].frontier.reset();
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

/** One non-empty fragment union in a duplicate class's sparse mask row. */
struct duplicateFragmentMaskEntry
{
    size_t fragment;
    const EdgeMaskAccumulator &mask;
};

/** Metadata stored separately from the level-wide accumulator buffer. */
struct duplicateFragmentMaskRow
{
    uint32_t fragment;
};

/**
 * @brief Read-only sparse fragment-mask row backed by one class level.
 *
 * Entries are sorted by fragment index. Iteration visits only fragments which
 * contain an occurrence, while fragmentCount retains the dense assembly-state
 * domain required by child enumeration and bound calculations.
 */
class duplicateFragmentMaskList
{
public:
    class iterator
    {
    public:
        using difference_type = ptrdiff_t;
        using value_type = duplicateFragmentMaskEntry;
        using iterator_category = forward_iterator_tag;

        iterator &operator++() noexcept
        {
            ++position;
            return *this;
        }

        iterator operator++(int) noexcept
        {
            iterator previous = *this;
            ++*this;
            return previous;
        }

        bool operator==(const iterator &other) const noexcept
        {
            return position == other.position;
        }

        duplicateFragmentMaskEntry operator*() const
        {
            return {
                rows[position].fragment,
                (*accumulators)[accumulatorOffset + position]
            };
        }

    private:
        friend class duplicateFragmentMaskList;

        iterator(
            span<const duplicateFragmentMaskRow> _rows,
            const EdgeMaskAccumulatorBuffer *_accumulators,
            size_t _accumulatorOffset,
            size_t _position
        ):
            rows(_rows),
            accumulators(_accumulators),
            accumulatorOffset(_accumulatorOffset),
            position(_position) {}

        span<const duplicateFragmentMaskRow> rows;
        const EdgeMaskAccumulatorBuffer *accumulators = nullptr;
        size_t accumulatorOffset = 0;
        size_t position = 0;
    };

    size_t fragmentCount = 0;

    bool empty() const noexcept {return rows.empty();}
    size_t size() const noexcept {return rows.size();}

    iterator begin() const
    {
        return iterator(rows, accumulators, accumulatorOffset, 0);
    }

    iterator end() const
    {
        return iterator(rows, accumulators, accumulatorOffset, rows.size());
    }

    duplicateFragmentMaskEntry operator[](size_t position) const
    {
        return *iterator(rows, accumulators, accumulatorOffset, position);
    }

    size_t maskCount(size_t fragment) const
    {
        const auto entry = lower_bound(
            rows.begin(),
            rows.end(),
            fragment,
            [](const duplicateFragmentMaskRow &candidate, size_t value)
            {
                return candidate.fragment < value;
            }
        );
        if (entry == rows.end() || entry->fragment != fragment) return 0;
        const size_t position = static_cast<size_t>(entry - rows.begin());
        return (*accumulators)[accumulatorOffset + position].count();
    }

private:
    template<typename> friend struct duplicateClassLevel;

    duplicateFragmentMaskList(
        size_t _fragmentCount,
        span<const duplicateFragmentMaskRow> _rows,
        const EdgeMaskAccumulatorBuffer &_accumulators,
        size_t _accumulatorOffset
    ):
        fragmentCount(_fragmentCount),
        rows(_rows),
        accumulators(&_accumulators),
        accumulatorOffset(_accumulatorOffset) {}

    span<const duplicateFragmentMaskRow> rows;
    const EdgeMaskAccumulatorBuffer *accumulators = nullptr;
    size_t accumulatorOffset = 0;
};

/** @brief One compact canonical-class bucket in an enumeration level. */
template<typename DuplicateSetType>
struct duplicateClassEntry
{
    int canonicalId;
    DuplicateSetType duplicates;

    explicit duplicateClassEntry(int _canonicalId): canonicalId(_canonicalId) {}

private:
    template<typename> friend struct duplicateClassLevel;

    static constexpr uint32_t noOccurrence =
        numeric_limits<uint32_t>::max();
    uint32_t firstOccurrence = noOccurrence;
    uint32_t lastOccurrence = noOccurrence;
    uint32_t occurrenceCount = 0;
    size_t maskOffset = 0;
    size_t maskCount = 0;
};

/**
 * @brief Compact duplicate classes with flat occurrences and CSR mask rows.
 *
 * Insertions append to one level-wide staging vector and link occurrences by
 * primitive indices. seal() sorts class metadata by canonical ID, moves every
 * class into one stable occurrence vector without changing its insertion
 * order, and constructs sparse fragment-union rows in flat accumulator storage.
 */
template<typename DuplicateSetType>
struct duplicateClassLevel
{
    using occurrence_type = typename DuplicateSetType::occurrence_type;
    using entry_type = duplicateClassEntry<DuplicateSetType>;

    class appender
    {
    public:
        void insert(occurrence_type occurrence)
        {
            level->insert(position, std::move(occurrence));
        }

    private:
        friend struct duplicateClassLevel;

        appender(duplicateClassLevel &_level, uint32_t _position):
            level(&_level), position(_position) {}

        duplicateClassLevel *level;
        uint32_t position;
    };

    vector<entry_type> classes;
    vector<uint8_t> aliveScratch;

    bool empty() const noexcept {return classes.empty();}
    size_t size() const noexcept {return classes.size();}
    size_t fragmentCount() const noexcept {return fragmentCountValue;}

    duplicateFragmentMaskList fragmentMasks(const entry_type &entry) const
    {
        const duplicateFragmentMaskRow *rowData = entry.maskCount == 0
            ? nullptr
            : maskRows.data() + entry.maskOffset;
        return duplicateFragmentMaskList(
            fragmentCountValue,
            span<const duplicateFragmentMaskRow>(rowData, entry.maskCount),
            maskAccumulators,
            entry.maskOffset
        );
    }

    /** Clear one reusable level without releasing any vector capacity. */
    void reset()
    {
        classes.clear();
        stagedOccurrences.clear();
        sealedOccurrences.clear();
        occurrenceOffsets.clear();
        maskRows.clear();
        maskOffsets.clear();
        fragmentPositions.clear();
        aliveScratch.clear();
        maskAccumulators.reset(0);
        duplicateSizeValue = 0;
        fragmentCountValue = 0;
        configured = false;
        fragmentCountPreset = false;
        sealed = false;
    }

    /** Reset and retain the assembly fragment domain for an empty level too. */
    void reset(size_t fragments)
    {
        reset();
        fragmentCountValue = fragments;
        fragmentCountPreset = true;
    }

    void seal()
    {
        if (sealed) return;
        if (classes.size() > 1)
        {
            sort(
                classes.begin(),
                classes.end(),
                [](const entry_type &left, const entry_type &right)
                {
                    return left.canonicalId < right.canonicalId;
                }
            );
        }

        compactOccurrences();
        buildFragmentMasks();
        stagedOccurrences.clear();
        sealed = true;
    }

    appender appendTo(uint32_t position)
    {
        if (sealed) throw logic_error("sealed duplicate-class level modified");
        if (position >= classes.size())
            throw logic_error("duplicate-class position is outside its level");
        return appender(*this, position);
    }

    void prepare(size_t duplicateSize, size_t fragments)
    {
        if (sealed) throw logic_error("sealed duplicate-class level modified");
        if (!configured)
        {
            if (fragmentCountPreset && fragmentCountValue != fragments)
            {
                throw logic_error("inconsistent duplicate-class fragment count");
            }
            duplicateSizeValue = duplicateSize;
            fragmentCountValue = fragments;
            configured = true;
            return;
        }
        if (
            duplicateSizeValue != duplicateSize ||
            fragmentCountValue != fragments
        )
        {
            throw logic_error("inconsistent duplicate-class level shape");
        }
    }

private:
    struct stagedOccurrence
    {
        occurrence_type occurrence;
        uint32_t next = entry_type::noOccurrence;

        explicit stagedOccurrence(occurrence_type _occurrence):
            occurrence(std::move(_occurrence)) {}
    };

    vector<stagedOccurrence> stagedOccurrences;
    vector<occurrence_type> sealedOccurrences;
    vector<size_t> occurrenceOffsets;
    vector<duplicateFragmentMaskRow> maskRows;
    vector<size_t> maskOffsets;
    vector<size_t> fragmentPositions;
    EdgeMaskAccumulatorBuffer maskAccumulators;
    size_t duplicateSizeValue = 0;
    size_t fragmentCountValue = 0;
    bool configured = false;
    bool fragmentCountPreset = false;
    bool sealed = false;

    void insert(uint32_t classPosition, occurrence_type occurrence)
    {
        if (sealed) throw logic_error("sealed duplicate-class level modified");
        if (classPosition >= classes.size())
            throw logic_error("duplicate-class position is outside its level");
        if (
            occurrence.fragment < 0 ||
            static_cast<size_t>(occurrence.fragment) >= fragmentCountValue
        )
        {
            throw logic_error("duplicate occurrence fragment is outside its level");
        }
        if (stagedOccurrences.size() >= entry_type::noOccurrence)
            throw length_error("duplicate occurrences exceed index capacity");

        entry_type &entry = classes[classPosition];
        if (entry.occurrenceCount == entry_type::noOccurrence)
            throw length_error("duplicate class exceeds occurrence capacity");
        const uint32_t occurrenceIndex = static_cast<uint32_t>(
            stagedOccurrences.size()
        );
        stagedOccurrences.emplace_back(std::move(occurrence));
        if (entry.lastOccurrence == entry_type::noOccurrence)
        {
            entry.firstOccurrence = occurrenceIndex;
        }
        else
        {
            stagedOccurrences[entry.lastOccurrence].next = occurrenceIndex;
        }
        entry.lastOccurrence = occurrenceIndex;
        ++entry.occurrenceCount;
    }

    void compactOccurrences()
    {
        sealedOccurrences.clear();
        sealedOccurrences.reserve(stagedOccurrences.size());
        occurrenceOffsets.clear();
        occurrenceOffsets.reserve(classes.size() + 1);

        for (entry_type &entry : classes)
        {
            occurrenceOffsets.push_back(sealedOccurrences.size());
            uint32_t occurrence = entry.firstOccurrence;
            uint32_t observedCount = 0;
            while (occurrence != entry_type::noOccurrence)
            {
                if (occurrence >= stagedOccurrences.size())
                    throw logic_error("invalid staged duplicate occurrence");
                stagedOccurrence &node = stagedOccurrences[occurrence];
                sealedOccurrences.push_back(std::move(node.occurrence));
                occurrence = node.next;
                ++observedCount;
            }
            if (observedCount != entry.occurrenceCount)
                throw logic_error("incomplete staged duplicate class");
        }
        occurrenceOffsets.push_back(sealedOccurrences.size());
    }

    void buildFragmentMasks()
    {
        constexpr size_t noPosition = numeric_limits<size_t>::max();
        maskRows.clear();
        maskOffsets.clear();
        maskOffsets.reserve(classes.size() + 1);
        fragmentPositions.assign(fragmentCountValue, noPosition);

        for (size_t classPosition = 0;
             classPosition < classes.size();
             ++classPosition)
        {
            maskOffsets.push_back(maskRows.size());
            const size_t occurrenceBegin = occurrenceOffsets[classPosition];
            const size_t occurrenceEnd = occurrenceOffsets[classPosition + 1];
            for (size_t occurrence = occurrenceBegin;
                 occurrence < occurrenceEnd;
                 ++occurrence)
            {
                const size_t fragment = static_cast<size_t>(
                    sealedOccurrences[occurrence].fragment
                );
                if (fragmentPositions[fragment] != noPosition) continue;
                fragmentPositions[fragment] = maskRows.size();
                maskRows.push_back({static_cast<uint32_t>(fragment)});
            }

            const size_t maskBegin = maskOffsets.back();
            sort(
                maskRows.begin() + maskBegin,
                maskRows.end(),
                [](const auto &left, const auto &right)
                {
                    return left.fragment < right.fragment;
                }
            );
            for (size_t row = maskBegin; row < maskRows.size(); ++row)
                fragmentPositions[maskRows[row].fragment] = noPosition;
        }
        maskOffsets.push_back(maskRows.size());

        maskAccumulators.reset(maskRows.size());
        for (size_t classPosition = 0;
             classPosition < classes.size();
             ++classPosition)
        {
            const size_t maskBegin = maskOffsets[classPosition];
            const size_t maskEnd = maskOffsets[classPosition + 1];
            for (size_t row = maskBegin; row < maskEnd; ++row)
                fragmentPositions[maskRows[row].fragment] = row;

            const size_t occurrenceBegin = occurrenceOffsets[classPosition];
            const size_t occurrenceEnd = occurrenceOffsets[classPosition + 1];
            for (size_t occurrence = occurrenceBegin;
                 occurrence < occurrenceEnd;
                 ++occurrence)
            {
                const occurrence_type &candidate = sealedOccurrences[occurrence];
                maskAccumulators[
                    fragmentPositions[static_cast<size_t>(candidate.fragment)]
                ].add(candidate.mask);
            }
            for (size_t row = maskBegin; row < maskEnd; ++row)
                fragmentPositions[maskRows[row].fragment] = noPosition;
        }

        for (size_t classPosition = 0;
             classPosition < classes.size();
             ++classPosition)
        {
            entry_type &entry = classes[classPosition];
            entry.maskOffset = maskOffsets[classPosition];
            entry.maskCount = maskOffsets[classPosition + 1] - entry.maskOffset;
            const size_t occurrenceBegin = occurrenceOffsets[classPosition];
            const size_t occurrenceCount =
                occurrenceOffsets[classPosition + 1] - occurrenceBegin;
            occurrence_type *occurrenceData = occurrenceCount == 0
                ? nullptr
                : sealedOccurrences.data() + occurrenceBegin;
            const bool spansFragments = entry.maskCount > 1;
            const bool valid = spansFragments || (
                entry.maskCount == 1 &&
                maskAccumulators[entry.maskOffset].count() >=
                    (duplicateSizeValue << 1)
            );
            entry.duplicates.bind(
                duplicateSizeValue,
                fragmentCountValue,
                span<occurrence_type>(occurrenceData, occurrenceCount),
                valid,
                spansFragments
            );
        }
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
    typename duplicateClassLevel<DuplicateSetType>::appender getOrCreate(
        duplicateClassLevel<DuplicateSetType> &level,
        int canonicalId,
        size_t duplicateSize,
        size_t fragmentCount
    )
    {
        if (generation == 0) [[unlikely]]
            throw logic_error("duplicate-class level was not started");
        level.prepare(duplicateSize, fragmentCount);
        if (
            !level.classes.empty() &&
            level.classes.back().canonicalId == canonicalId
        ) [[likely]]
        {
            return level.appendTo(
                static_cast<uint32_t>(level.classes.size() - 1)
            );
        }
        const size_t id = static_cast<size_t>(canonicalId);
        if (canonicalId >= 0 && id < slots.size()) [[likely]]
        {
            const slot &entry = slots[id];
            if (entry.generation == generation) [[likely]]
            {
                return level.appendTo(entry.position);
            }
        }
        return create(
            level,
            canonicalId
        );
    }

private:
    template<typename DuplicateSetType>
    [[gnu::noinline]]
    typename duplicateClassLevel<DuplicateSetType>::appender create(
        duplicateClassLevel<DuplicateSetType> &level,
        int canonicalId
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
        level.classes.emplace_back(canonicalId);
        return level.appendTo(entry.position);
    }
};

/**
 * @brief Generate the next set of duplicates from the duplicate d
 * 
 * @param dag Runtime DAG used to enumerate child masks
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
    const vector<dagLevel> &dag,
    potentialDuplicate &d,
    dagDuplicateClassLevel &stmap,
    duplicateClassIndexWorkspace &classIndex,
    const EdgeMask &fragment,
    size_t size,
    int ordinal,
    size_t frags
)
{
    bool overweight = 0;
    const dagLevel &parentLevel = dag[size - 1];
    const dagNode &parent = parentLevel.nodes[d.idx];
    if (parent.transitionCount == 0) return overweight;
    const dagTransition *transition =
        parentLevel.transitions.data() + parent.transitionOffset;
    const dagTransition *const transitionEnd =
        transition + parent.transitionCount;
    const dagNode *const childNodes = dag[size].nodes.data();
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
 * @param dag Runtime DAG used to enumerate child masks
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
    const vector<dagLevel> &dag,
    dagDuplicateSet &ds,
    dagDuplicateClassLevel &stmap,
    duplicateClassIndexWorkspace &classIndex,
    span<EdgeMaskAccumulator> takenMasks,
    const vector<assemblyFragment> &stateFragments,
    int ordinal,
    bool &overweight,
    bool last,
    vector<uint8_t> &aliveScratch
)
    {
        bool output = 0;
        const bool allAlive = ds.occurrencesSpanMultipleFragments();
        if (searchShouldStop()) return output;

        auto generateFromDuplicate = [&](potentialDuplicate &duplicate)
        {
            const int frag = duplicate.fragment;
            takenMasks[frag].add(duplicate.mask);
            ds.dead = 0;
            if (!last)
            {
                overweight |= dagGenerate(
                    dag,
                    duplicate,
                    stmap,
                    classIndex,
                    stateFragments[frag].mask,
                    ds.size,
                    ordinal,
                    ds.fragmentCount
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

        aliveScratch.assign(ds.list.size(), 0);
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
                        aliveScratch[i] = 1;
                        aliveScratch[j] = 1;
                    }
                }
                else
                {
                    aliveScratch[i] = 1;
                    aliveScratch[j] = 1;
                }
            }
            if (aliveScratch[i])
            {
                if (!generateFromDuplicate(ds.list[i])) return output;
            }
        }
        return output;
    }
