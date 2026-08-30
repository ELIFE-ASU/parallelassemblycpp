/**
 * @brief Enumerate all subgraphs during the initial phase of the pathway algorithm. See Seet et al section 4.3 Duplicate Enumeration
 *
 * @param _target The initial assembly state
 * @param stmapVector The matchings found
 * @return true if any matchings found
 * @return false if no matchings found
 */
bool initialRecursiveEnumeration(
    assemblyState &_target,
    vector<initialDuplicateClassLevel> &stmapVector,
    duplicateClassIndexWorkspace &classIndex,
    vector<dagLevel> &dag
)
{
    vector<initialDagLevel> tempDag(2);
    vector<assemblyFragment> &fragments = _target.fragments;
    const initialIncidentEdgeIndex incidentEdges;
    bool alive = 0;
    size_t currSize = 1;
    vector<initialPotentialDuplicate> prevML;
    vector<int> rootNodeIndices(univEdgeList.size(), -1);
    size_t retainedStateCount = 0;

    // Retain the one-edge DAG states first so they count toward ENUM_MAX too.
    for (size_t i = 0; i < fragments.size(); i++)
    {
        if (searchShouldStop()) return false;
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStopPeriodically()) return false;
            if (fragments[i].mask[j] != 0)
            {
                EdgeMask b = 0; b.set(j);
                const initialDagInsertion insertion =
                    tryRetainInitialDagMask(
                        tempDag[0],
                        b,
                        retainedStateCount
                    );
                if (insertion.result == initialDagInsertionResult::limitReached)
                    return false;
                rootNodeIndices[j] = insertion.index;
            }
        }
    }

    // Generate the first multi-edge frontier only after all base states exist.
    for (size_t i = 0; i < fragments.size(); i++)
    {
        if (searchShouldStop()) return false;
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStopPeriodically()) return false;
            if (fragments[i].mask[j] == 0) continue;

            if (rootNodeIndices[j] < 0)
                throw logic_error("initial DAG root was not retained");
            initialPotentialDuplicate m(
                j,
                fragments[i].mask,
                incidentEdges,
                i,
                rootNodeIndices[j]
            );
            if (!m.generateDAG(
                prevML,
                retainedStateCount,
                tempDag,
                fragments[i].mask,
                incidentEdges
            ))
            {
                return false;
            }
        }
    }

    bool active = true;
    while (active)
    {
        if (searchShouldStop()) return false;
        stmapVector.emplace_back();
        initialDuplicateClassLevel &stmap = stmapVector.back();
        classIndex.beginLevel();
        active = 0;
        vector<initialPotentialDuplicate> currML;
        for (size_t i = 0; i < prevML.size(); i++)
        {
            if (searchShouldStop()) return false;
            initialPotentialDuplicate &m = prevML[i];

            int s = canonise(m.mask);
            if (searchShouldStop()) return false;
            classIndex.getOrCreate(
                stmap,
                s,
                currSize + 1,
                fragments.size()
            ).insert(std::move(m));
        }
        stmap.seal();
        tempDag.resize(tempDag.size() + 1);
        for (auto &entry : stmap.classes)
        {
            if (searchShouldStop()) return false;
            initialDuplicateSet &ss = entry.duplicates;
            if (ss.isValid())
            {
                if (searchShouldStop()) return false;
                active = 1;
                alive |= ss.dagPopulator(
                    currML,
                    retainedStateCount,
                    tempDag,
                    fragments,
                    incidentEdges,
                    stmap.aliveScratch
                );
                if (enumerationLimitReached)
                {
                    return false;
                }
                if (searchShouldStop()) return false;
            }
        }
        currSize++;
        prevML = std::move(currML);
    }
    if (searchShouldStop()) return false;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    setSearchTelemetryPhase(SearchTelemetryPhase::dagConversion);
#endif
    convertDag(tempDag, dag);
    if (searchShouldStop()) return false;
    return alive;
}

/**
 * @brief Enumerate all subgraphs during subsequent phases of the pathway algorithm using the DAG to speed things up.  See seet et al section 4.3 Duplicate Enumeration
 *
 * @param _target Target assembly state
 * @param stmapVector List of generated duplicate pairs
 * @return Maximum duplicate size reached
 */
int dagRecursiveEnumeration(
    const vector<dagLevel> &dag,
    assemblyState &_target,
    dagAssemblySearchFrame &frame,
    duplicateClassIndexWorkspace &classIndex
)
{
    if (searchShouldStop()) return 0;
    int ordinal = std::numeric_limits<int>::max();
    if (_target.fragments.front().canonicalId >= 0)
        ordinal = _target.fragments.front().canonicalId;
    vector<assemblyFragment> &fragments = _target.fragments;
    size_t currSize = 1;
    
    dagDuplicateClassLevel &firstLevel =
        frame.appendDuplicateLevel(fragments.size());
    classIndex.beginLevel();
    for (size_t i = 0; i < fragments.size(); i++)
    {
        if (searchShouldStop()) return 0;
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStopPeriodically()) return 0;
            if (fragments[i].mask[j] != 0)
            {
                EdgeMask b = 0; b.set(j);
                potentialDuplicate m(std::move(b), i, j);
                dagGenerate(
                    dag,
                    m,
                    firstLevel,
                    classIndex,
                    fragments[i].mask,
                    currSize,
                    ordinal,
                    fragments.size()
                );
                if (searchShouldStop()) return 0;
            }
        }
    }
    firstLevel.seal();
    bool active = 1, overweight = 0, last = 0;
    while (active)
    {
        if (searchShouldStop()) return 0;
        span<EdgeMaskAccumulator> targetMask = frame.targetMasks.appendRow();
        active = 0;
        dagDuplicateClassLevel &nextLevel =
            frame.appendDuplicateLevel(fragments.size());
        classIndex.beginLevel();
        dagDuplicateClassLevel &stmap =
            frame.duplicateLevels[frame.duplicateLevelCount - 2];
        for (auto &entry : stmap.classes)
        {
            if (searchShouldStop()) return 0;
            dagDuplicateSet &ss = entry.duplicates;
            if (ss.isValid())
            {
                if (searchShouldStop()) return 0;
                active |= dagDuplicateGenerator(
                    dag,
                    ss,
                    nextLevel,
                    classIndex,
                    targetMask,
                    fragments,
                    ordinal,
                    overweight,
                    last,
                    stmap.aliveScratch
                );
                if (searchShouldStop()) return 0;
            }
        }
        nextLevel.seal();
        if (overweight) last = 1;
        currSize++;
    }
    if (
        frame.duplicateLevels[frame.duplicateLevelCount - 1].empty()
    ) --frame.duplicateLevelCount;
    for (size_t i = 0; i < fragments.size(); i++)
    {
        if (searchShouldStop()) return 0;
        assemblyFragment &fragment = fragments[i];
        fragment.retainEdges(frame.targetMasks.row(0)[i]);
    }
    return currSize;
}

/**
 * @brief Record a stable edge order for uniformly labelled path molecules.
 *
 * On such a graph, a connected fragment's canonical class is determined by
 * its edge count. The order lets matching-pair filtering describe residual
 * components without running union-find first.
 */
void configureHomogeneousPathEdgePositions(vector<int> &edgePositions)
{
    edgePositions.clear();
    const size_t atomCount = targetMolecule.mg.size();
    const size_t edgeCount = univEdgeList.size();
    if (edgeCount < 2 || atomCount != edgeCount + 1) return;

    const string &atomType = targetMolecule.mg.front().type;
    for (const atom &candidate : targetMolecule.mg)
    {
        if (candidate.type != atomType) return;
    }

    vector<vector<pair<int, int>>> adjacency(atomCount);
    short bondType = 0;
    bool foundBondType = false;
    for (size_t edge = 0; edge < edgeCount; edge++)
    {
        const edgeL &entry = univEdgeList[edge];
        const short candidateBondType = targetMolecule.btypeS(entry.a, entry.c);
        if (!foundBondType)
        {
            bondType = candidateBondType;
            foundBondType = true;
        }
        else if (candidateBondType != bondType)
            return;
        adjacency[entry.a].emplace_back(entry.b, static_cast<int>(edge));
        adjacency[entry.b].emplace_back(entry.a, static_cast<int>(edge));
    }

    size_t endpoint = atomCount;
    size_t endpointCount = 0;
    for (size_t atomIndex = 0; atomIndex < atomCount; atomIndex++)
    {
        const size_t degree = adjacency[atomIndex].size();
        if (degree == 1)
        {
            endpoint = atomIndex;
            ++endpointCount;
        }
        else if (degree != 2)
        {
            return;
        }
    }
    if (endpointCount != 2) return;

    vector<int> candidatePositions(edgeCount, -1);
    size_t previous = atomCount;
    size_t current = endpoint;
    for (size_t position = 0; position < edgeCount; position++)
    {
        const auto &neighbours = adjacency[current];
        const auto next = neighbours.front().first != static_cast<int>(previous)
            ? neighbours.front()
            : neighbours.back();
        if (candidatePositions[next.second] != -1) return;
        candidatePositions[next.second] = static_cast<int>(position);
        previous = current;
        current = static_cast<size_t>(next.first);
    }
    if (adjacency[current].size() != 1) return;
    edgePositions = std::move(candidatePositions);
}


/**
 * @brief Compute the targeted and unrestricted post-fragment bounds together
 * in one pass over the child fragments.
 * @param target The target assembly state
 * @param matchMask The bitset of all graphs isomorphic to the matching
 * @param maxFragMask The bitset of all duplicatable subgraphs with the same bitset count as the matching
 * @param boundTotals Reusable scratch for each unrestricted duplicate size
 * @return int the largest remaining duplicate-bond estimate
 */
template<typename MatchMask, typename MaxFragmentMask>
int postFragmentationCutoff(
    const assemblyState &target,
    const MatchMask &matchMask,
    const MaxFragmentMask &maxFragMask,
    vi &boundTotals
)
{
    const int maxFragSize = target.fragments[0].edgeCount;
    const bool useTargetedBounds = target.fragments.size() >= 2;
    int matchDB = 0;
    int maxFragDB = 0;
    boundTotals.assign(max(maxFragSize - 2, 0), 0);

    for (size_t i = 0; i < target.fragments.size(); i++)
    {
        const assemblyFragment &fragment = target.fragments[i];
        const int edgeCount = fragment.edgeCount;
        if (useTargetedBounds)
        {
            const int matchingEdges = i == 0
                ? maxFragSize
                : static_cast<int>(
                    matchMask.intersectionCount(fragment.mask)
                );
            const int maximalEdges = i == 0
                ? maxFragSize
                : static_cast<int>(
                    maxFragMask.intersectionCount(fragment.mask)
                );
            matchDB += assemblyState::fixedSizeDupBondsForFragment(
                edgeCount,
                maxFragSize,
                matchingEdges
            );
            maxFragDB += assemblyState::fixedSizeDupBondsForFragment(
                edgeCount,
                maxFragSize,
                maximalEdges
            );
        }

        if (!boundTotals.empty()) boundTotals[0] += edgeCount / 2;
        for (int duplicateSize = 3;
             duplicateSize < maxFragSize;
             duplicateSize++)
        {
            boundTotals[duplicateSize - 2] +=
                assemblyState::unrestrictedDupBondsForFragment(
                    edgeCount,
                    duplicateSize
                );
        }
    }

    int result = 0;
    if (useTargetedBounds)
    {
        matchDB -= ceilLog2(maxFragSize);
        maxFragDB -= ceilLog2(maxFragSize) + 1;
        result = max(matchDB, maxFragDB);
    }
    if (!boundTotals.empty())
    {
        boundTotals[0]--;
        result = max(result, boundTotals[0]);
        for (int duplicateSize = 3;
             duplicateSize < maxFragSize;
             duplicateSize++)
        {
            const size_t index = duplicateSize - 2;
            boundTotals[index] -= ceilLog2(duplicateSize);
            result = max(result, boundTotals[index]);
        }
    }
    return result;
}

void buildUnrestrictedDupBondTotals(
    const assemblyState &target,
    int maxDuplicateSize,
    vi &totals
)
{
    totals.assign(max(maxDuplicateSize - 1, 0), 0);
    for (const assemblyFragment &fragment : target.fragments)
    {
        const int edgeCount = fragment.edgeCount;
        if (!totals.empty()) totals[0] += edgeCount / 2;
        for (int duplicateSize = 3;
             duplicateSize <= maxDuplicateSize;
             duplicateSize++)
        {
            totals[duplicateSize - 2] +=
                assemblyState::unrestrictedDupBondsForFragment(
                    edgeCount,
                    duplicateSize
                );
        }
    }
}

int pairSpecificGenericBound(
    const assemblyState &target,
    const validMatchings &matching,
    const vi &parentTotals
)
{
    // Evaluate the generic bound on a virtual child whose selected copies have
    // been removed but whose residual parents remain unsplit. Since
    // n - ceil(n / k) is superadditive, later component splitting can only
    // lower this duplicate-bond estimate, so it is safe before union-find.
    const int selectedSize = matching.maxFragSize;
    int total = parentTotals[0] + selectedSize / 2;
    if (matching.frag1 == matching.frag2)
    {
        const int parentEdges =
            target.fragments[matching.frag1].edgeCount;
        total += (parentEdges - 2 * selectedSize) / 2 - parentEdges / 2;
    }
    else
    {
        const int firstParentEdges =
            target.fragments[matching.frag1].edgeCount;
        const int secondParentEdges =
            target.fragments[matching.frag2].edgeCount;
        total += (firstParentEdges - selectedSize) / 2 - firstParentEdges / 2;
        total += (secondParentEdges - selectedSize) / 2 - secondParentEdges / 2;
    }
    int result = total - 1;

    for (int duplicateSize = 3;
         duplicateSize < selectedSize;
         duplicateSize++)
    {
        total = parentTotals[duplicateSize - 2] +
            assemblyState::unrestrictedDupBondsForFragment(
                selectedSize,
                duplicateSize
            );
        if (matching.frag1 == matching.frag2)
        {
            const int parentEdges =
                target.fragments[matching.frag1].edgeCount;
            total += assemblyState::unrestrictedDupBondsForFragment(
                parentEdges - 2 * selectedSize,
                duplicateSize
            ) - assemblyState::unrestrictedDupBondsForFragment(
                parentEdges,
                duplicateSize
            );
        }
        else
        {
            const int firstParentEdges =
                target.fragments[matching.frag1].edgeCount;
            const int secondParentEdges =
                target.fragments[matching.frag2].edgeCount;
            total += assemblyState::unrestrictedDupBondsForFragment(
                firstParentEdges - selectedSize,
                duplicateSize
            ) - assemblyState::unrestrictedDupBondsForFragment(
                firstParentEdges,
                duplicateSize
            );
            total += assemblyState::unrestrictedDupBondsForFragment(
                secondParentEdges - selectedSize,
                duplicateSize
            ) - assemblyState::unrestrictedDupBondsForFragment(
                secondParentEdges,
                duplicateSize
            );
        }
        result = max(result, total - ceilLog2(duplicateSize));
    }
    return result;
}

[[gnu::noinline]] int pairSpecificGenericBound(
    const assemblyState &target,
    int selectedSize,
    int firstFragment,
    int secondFragment,
    const vi &parentTotals
)
{
    // Evaluate the generic bound on a virtual child whose selected copies have
    // been removed but whose residual parents remain unsplit. Since
    // n - ceil(n / k) is superadditive, later component splitting can only
    // lower this duplicate-bond estimate, so it is safe before union-find.
    int total = parentTotals[0] + selectedSize / 2;
    if (firstFragment == secondFragment)
    {
        const int parentEdges =
            target.fragments[firstFragment].edgeCount;
        total += (parentEdges - 2 * selectedSize) / 2 - parentEdges / 2;
    }
    else
    {
        const int firstParentEdges =
            target.fragments[firstFragment].edgeCount;
        const int secondParentEdges =
            target.fragments[secondFragment].edgeCount;
        total += (firstParentEdges - selectedSize) / 2 - firstParentEdges / 2;
        total += (secondParentEdges - selectedSize) / 2 - secondParentEdges / 2;
    }
    int result = total - 1;

    for (int duplicateSize = 3;
         duplicateSize < selectedSize;
         duplicateSize++)
    {
        total = parentTotals[duplicateSize - 2] +
            assemblyState::unrestrictedDupBondsForFragment(
                selectedSize,
                duplicateSize
            );
        if (firstFragment == secondFragment)
        {
            const int parentEdges =
                target.fragments[firstFragment].edgeCount;
            total += assemblyState::unrestrictedDupBondsForFragment(
                parentEdges - 2 * selectedSize,
                duplicateSize
            ) - assemblyState::unrestrictedDupBondsForFragment(
                parentEdges,
                duplicateSize
            );
        }
        else
        {
            const int firstParentEdges =
                target.fragments[firstFragment].edgeCount;
            const int secondParentEdges =
                target.fragments[secondFragment].edgeCount;
            total += assemblyState::unrestrictedDupBondsForFragment(
                firstParentEdges - selectedSize,
                duplicateSize
            ) - assemblyState::unrestrictedDupBondsForFragment(
                firstParentEdges,
                duplicateSize
            );
            total += assemblyState::unrestrictedDupBondsForFragment(
                secondParentEdges - selectedSize,
                duplicateSize
            ) - assemblyState::unrestrictedDupBondsForFragment(
                secondParentEdges,
                duplicateSize
            );
        }
        result = max(result, total - ceilLog2(duplicateSize));
    }
    return result;
}

template<typename DuplicateMasks>
[[gnu::noinline]] bool shouldVisitFragmentPairBlocks(
    assemblyState &target,
    const dagDuplicateSet &duplicates,
    const DuplicateMasks &duplicateMasks,
    int targetedBound,
    int pairBoundLimit,
    int &matchingClassBound
)
{
    // In small classes only a single fragment run has enough reuse to cover
    // the block-dispatch cost. The size-three generic bound has no inner
    // duplicate-size loop, so keep that rule for the minimum 16-occurrence
    // block as well. The endpoint check rejects these common multi-run cases
    // before the full sortedness scan.
    if (
        (
            duplicates.list.size() < 16 ||
            (
                duplicates.list.size() == 16 &&
                duplicates.size == 3
            )
        ) &&
        duplicates.list.front().fragment !=
            duplicates.list.back().fragment
    ) return false;

    if (targetedBound - 1 > pairBoundLimit) return false;
    if (targetedBound > pairBoundLimit)
    {
        matchingClassBound = target.maxDupBonds(
            duplicates.size,
            duplicateMasks
        );
        if (matchingClassBound > pairBoundLimit) return false;
    }
    return duplicates.hasDenseFragmentRuns();
}

struct homogeneousPathResidualKey
{
    // At most two removed parents and four residual path components.
    array<int, 12> values{};
    unsigned char used = 0;

    bool operator==(const homogeneousPathResidualKey &other) const
    {
        return
            used == other.used &&
            equal(values.begin(), values.begin() + used, other.values.begin());
    }
};

struct homogeneousPathResidualKeyHash
{
    size_t operator()(const homogeneousPathResidualKey &key) const
    {
        size_t result = key.used;
        for (size_t i = 0; i < key.used; i++)
        {
            result ^= static_cast<size_t>(key.values[i]) + 0x9e3779b9 +
                (result << 6) + (result >> 2);
        }
        return result;
    }
};

/**
 * @brief Quotient matching pairs on a uniformly labelled path by child state.
 *
 * Every connected fragment is an interval in the molecule-wide edge order.
 * Its canonical class is therefore determined solely by its edge count. A
 * pair's exact child is represented by the small multiset delta obtained by
 * removing its parent path(s) and adding the residual intervals.
 */
template<typename DuplicateSet>
struct homogeneousPathEquivalentMatchings
{
    const assemblyState &input;
    const DuplicateSet &duplicates;
    const vector<int> &edgePositions;
    bool enabled = false;
    vector<int> fragmentStarts;
    vector<int> occurrenceStarts;
    unordered_set<
        homogeneousPathResidualKey,
        homogeneousPathResidualKeyHash
    > seen;

    homogeneousPathEquivalentMatchings(
        const assemblyState &_input,
        const DuplicateSet &_duplicates,
        const vector<int> &_edgePositions
    ):
        input(_input),
        duplicates(_duplicates),
        edgePositions(_edgePositions)
    {
        constexpr uint64_t minimumValidPairs = 16;
        if (
            edgePositions.empty() ||
            duplicates.list.size() < 2 ||
            duplicates.size > static_cast<size_t>(numeric_limits<int>::max())
        ) return;

        fragmentStarts.resize(input.fragments.size());
        for (size_t fragment = 0; fragment < input.fragments.size(); fragment++)
        {
            if (!intervalStart(
                input.fragments[fragment].mask,
                input.fragments[fragment].edgeCount,
                fragmentStarts[fragment]
            )) return;
        }

        occurrenceStarts.resize(duplicates.list.size());
        vector<vector<int>> startsByFragment(input.fragments.size());
        for (size_t occurrence = 0;
             occurrence < duplicates.list.size();
             occurrence++)
        {
            const auto &candidate = duplicates.list[occurrence];
            if (
                candidate.fragment < 0 ||
                static_cast<size_t>(candidate.fragment) >=
                    input.fragments.size() ||
                !input.fragments[candidate.fragment].mask.contains(
                    candidate.mask
                ) ||
                !intervalStart(
                    candidate.mask,
                    static_cast<int>(duplicates.size),
                    occurrenceStarts[occurrence]
                )
            ) return;
            startsByFragment[candidate.fragment].push_back(
                occurrenceStarts[occurrence]
            );
        }

        uint64_t validPairs = 0;
        uint64_t priorOccurrences = 0;
        const int duplicateSize = static_cast<int>(duplicates.size);
        for (vector<int> &starts : startsByFragment)
        {
            if (starts.empty()) continue;
            validPairs += min<uint64_t>(
                minimumValidPairs - min(validPairs, minimumValidPairs),
                priorOccurrences * starts.size()
            );
            if (validPairs >= minimumValidPairs) break;
            priorOccurrences += starts.size();

            sort(starts.begin(), starts.end());
            for (size_t first = 0; first + 1 < starts.size(); first++)
            {
                const auto second = lower_bound(
                    starts.begin() + first + 1,
                    starts.end(),
                    starts[first] + duplicateSize
                );
                validPairs += static_cast<uint64_t>(starts.end() - second);
                if (validPairs >= minimumValidPairs) break;
            }
            if (validPairs >= minimumValidPairs) break;
        }
        if (validPairs < minimumValidPairs) return;
        seen.reserve(256);
        enabled = true;
    }

    bool skip(
        const validMatchings &matching,
        size_t firstOccurrence,
        size_t secondOccurrence
    )
    {
        if (!enabled) return false;
        array<pair<int, int>, 6> changes;
        size_t changeCount = 0;
        const int duplicateSize = static_cast<int>(duplicates.size);
        const auto addChange = [&](int edgeCount, int delta)
        {
            if (edgeCount >= 2)
                changes[changeCount++] = {edgeCount, delta};
        };

        if (matching.frag1 == matching.frag2)
        {
            const int fragment = matching.frag1;
            int firstStart = occurrenceStarts[firstOccurrence] -
                fragmentStarts[fragment];
            int secondStart = occurrenceStarts[secondOccurrence] -
                fragmentStarts[fragment];
            if (secondStart < firstStart) swap(firstStart, secondStart);
            if (secondStart < firstStart + duplicateSize) return false;

            const int parentEdges = input.fragments[fragment].edgeCount;
            addChange(parentEdges, -1);
            addChange(firstStart, 1);
            addChange(secondStart - firstStart - duplicateSize, 1);
            addChange(parentEdges - secondStart - duplicateSize, 1);
        }
        else
        {
            const int fragments[2] = {matching.frag1, matching.frag2};
            const size_t occurrences[2] = {
                firstOccurrence,
                secondOccurrence
            };
            for (size_t selected = 0; selected < 2; selected++)
            {
                const int fragment = fragments[selected];
                const int start = occurrenceStarts[occurrences[selected]] -
                    fragmentStarts[fragment];
                const int parentEdges = input.fragments[fragment].edgeCount;
                addChange(parentEdges, -1);
                addChange(start, 1);
                addChange(parentEdges - start - duplicateSize, 1);
            }
        }

        sort(changes.begin(), changes.begin() + changeCount);
        homogeneousPathResidualKey key;
        for (size_t index = 0; index < changeCount;)
        {
            const int edgeCount = changes[index].first;
            int delta = 0;
            do
            {
                delta += changes[index].second;
                ++index;
            }
            while (
                index < changeCount &&
                changes[index].first == edgeCount
            );
            if (delta == 0) continue;
            key.values[key.used++] = edgeCount;
            key.values[key.used++] = delta;
        }
        return !seen.insert(key).second;
    }

private:
    bool intervalStart(
        const EdgeMask &mask,
        int expectedEdgeCount,
        int &start
    ) const
    {
        if (expectedEdgeCount < 1) return false;
        int minimum = numeric_limits<int>::max();
        int maximum = -1;
        int visited = 0;
        for (size_t edge = mask.findFirst();
             edge < edgePositions.size();
             edge = mask.findNext(edge))
        {
            const int position = edgePositions[edge];
            if (position < 0) return false;
            minimum = min(minimum, position);
            maximum = max(maximum, position);
            ++visited;
        }
        if (
            visited != expectedEdgeCount ||
            maximum - minimum + 1 != expectedEdgeCount
        ) return false;
        start = minimum;
        return true;
    }
};

enum class matchingEquivalenceMode
{
    none,
    homogeneousPath
};

enum class fragmentPairBlockMode
{
    automatic,
    disabled,
    enabled
};

// Pair pruning pays for itself once residual splitting is no longer dominated
// by the fixed cost of the bound lookup.
constexpr size_t pairBoundMinimumMoleculeEdges = 27;

template<
    matchingEquivalenceMode equivalenceMode,
    bool trackPath,
    bool allowDepthTwoDonation = false,
    bool useSharedStates = false,
    bool enableFragmentPairBlocks = false
>
void dagRecursiveAssemblyWithWorkspaceImpl(
    const vector<dagLevel> &dag,
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
);

template<bool trackPath>
void recordImprovedAssemblyIndex(
    assemblyState &input,
    int &AI,
    const assemblySearchStorage &searchStorage
)
{
    if (sharedAssemblyIndex != nullptr)
        AI = min(AI, sharedAssemblyIndex->load(std::memory_order_relaxed));
    const int candidate = input.AI();
    if (candidate >= AI) return;

    AI = candidate;
    if (sharedAssemblyIndex != nullptr)
    {
        int observed = sharedAssemblyIndex->load(std::memory_order_relaxed);
        while (
            candidate < observed &&
            !sharedAssemblyIndex->compare_exchange_weak(
                observed,
                candidate,
                std::memory_order_relaxed
            )
        ) {}
        AI = min(AI, observed);
        if (candidate > AI) return;
    }
    if constexpr (trackPath)
        searchStorage.pathway->best = searchStorage.pathway->current;
    const unsigned long long time = elapsedClockTicks();
    if (!suppressSearchOutput)
    {
#ifdef ASSEMBLYCPP_LIBRARY_BUILD
        if (verbose)
#endif
        cout << "Best assembly index: " << AI << " (" << time
             << " clock ticks)\n";
    }
    if (writeIntermediateMAs) intermediateMAs.emplace_back(time, AI);
}

template<
    matchingEquivalenceMode equivalenceMode,
    bool trackPath,
    bool allowDepthTwoDonation = false,
    bool useSharedStates = false,
    fragmentPairBlockMode pairBlockMode = fragmentPairBlockMode::automatic
>
bool continueCanonicalAssemblySearchWithWorkspace(
    const vector<dagLevel> &dag,
    assemblyState &candidate,
    span<const int> candidateKey,
    int sumDupBonds,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage,
    validMatchings *matching = nullptr
)
{
    if (searchShouldStop()) return false;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.assemblyCacheLookups;
#endif

    const assemblyTranspositionTable::result result = [&]
    {
        if constexpr (useSharedStates)
            return searchStorage.considerShared(candidateKey, sumDupBonds);
        else
            return searchStorage.states.consider(candidateKey, sumDupBonds);
    }();

#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
    {
        const bool inserted =
            result == assemblyTranspositionTable::result::inserted;
        if (inserted) ++searchTelemetry.counters.assemblyCacheMisses;
        else ++searchTelemetry.counters.assemblyCacheHits;
    }
#endif

    if (result == assemblyTranspositionTable::result::dominated)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.assemblyCachePrunedHits;
#endif
        return true;
    }
    if (result == assemblyTranspositionTable::result::improved)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.assemblyCacheUpdatedHits;
#endif
    }

    if constexpr (trackPath)
    {
        if (matching == nullptr)
            throw logic_error("path-tracking search is missing its matching");
        searchStorage.pathway->current.push_back(
            assemblyPathStep{matching->first, matching->second}
        );
    }
    if constexpr (
        equivalenceMode == matchingEquivalenceMode::none && !trackPath
    )
    {
        if constexpr (pairBlockMode == fragmentPairBlockMode::enabled)
        {
            dagRecursiveAssemblyWithWorkspaceImpl<
                equivalenceMode,
                trackPath,
                allowDepthTwoDonation,
                useSharedStates,
                true
            >(
                dag,
                candidate,
                AI,
                fragmentationWorkspace,
                searchStorage
            );
        }
        else if constexpr (pairBlockMode == fragmentPairBlockMode::disabled)
        {
            dagRecursiveAssemblyWithWorkspaceImpl<
                equivalenceMode,
                trackPath,
                allowDepthTwoDonation,
                useSharedStates,
                false
            >(
                dag,
                candidate,
                AI,
                fragmentationWorkspace,
                searchStorage
            );
        }
        else if (
            fragmentationWorkspace.edgeCount >=
                pairBoundMinimumMoleculeEdges &&
            fragmentationWorkspace.edgeCount <= EdgeMask::wordBits
        )
        {
            dagRecursiveAssemblyWithWorkspaceImpl<
                equivalenceMode,
                trackPath,
                allowDepthTwoDonation,
                useSharedStates,
                true
            >(
                dag,
                candidate,
                AI,
                fragmentationWorkspace,
                searchStorage
            );
        }
        else
        {
            dagRecursiveAssemblyWithWorkspaceImpl<
                equivalenceMode,
                trackPath,
                allowDepthTwoDonation,
                useSharedStates,
                false
            >(
                dag,
                candidate,
                AI,
                fragmentationWorkspace,
                searchStorage
            );
        }
    }
    else
    {
        dagRecursiveAssemblyWithWorkspaceImpl<
            equivalenceMode,
            trackPath,
            allowDepthTwoDonation,
            useSharedStates,
            false
        >(
            dag,
            candidate,
            AI,
            fragmentationWorkspace,
            searchStorage
        );
    }
    if constexpr (trackPath) searchStorage.pathway->current.pop_back();
    return !searchShouldStop();
}

template<
    matchingEquivalenceMode equivalenceMode,
    bool trackPath,
    bool allowDepthTwoDonation = false,
    bool useSharedStates = false,
    fragmentPairBlockMode pairBlockMode = fragmentPairBlockMode::automatic
>
bool continueAssemblySearchWithWorkspace(
    const vector<dagLevel> &dag,
    assemblyState &candidate,
    validMatchings &matching,
    span<const int> candidateKey,
    int sumDupBonds,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    return continueCanonicalAssemblySearchWithWorkspace<
        equivalenceMode,
        trackPath,
        allowDepthTwoDonation,
        useSharedStates,
        pairBlockMode
    >(
        dag,
        candidate,
        candidateKey,
        sumDupBonds,
        AI,
        fragmentationWorkspace,
        searchStorage,
        &matching
    );
}

/**
 * @brief The recursive function that enumerates duplicates and generates assembly states on all but the first pass
 * of the assembly algorithm
 * 
 * @param input The input assembly state
 * @param AI The global minimum assembly index found
 * @param fragmentationWorkspace Buffers reused across the search
 */
template<
    matchingEquivalenceMode equivalenceMode,
    bool trackPath,
    bool allowDepthTwoDonation,
    bool useSharedStates,
    bool enableFragmentPairBlocks
>
void dagRecursiveAssemblyWithWorkspaceImpl(
    const vector<dagLevel> &dag,
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    const bool usePairBound =
        fragmentationWorkspace.edgeCount >= pairBoundMinimumMoleculeEdges;
    recordImprovedAssemblyIndex<trackPath>(input, AI, searchStorage);
    if (searchShouldStop()) return;

    dagAssemblySearchFrameScope frameScope(
        searchStorage,
        input.fragments.size()
    );
    dagAssemblySearchFrame &frame = frameScope.frame;
    int maxFragSize = dagRecursiveEnumeration(
        dag,
        input,
        frame,
        searchStorage.duplicateClassIndex
    );

    if (searchShouldStop()) return;

    /// Find the fragment-size-specific AI lower bounds
    vi &fragSizeListMax = frame.fragSizeListMax;
    input.maxDupBondsPrefix(fragSizeListMax, maxFragSize, frame.targetMasks);
    if (searchShouldStop()) return;

    /// Begin iterating through the enumerated duplicatable fragments
    assemblyState &candidate = frame.candidate;
    for (size_t levelIndex = frame.duplicateLevelCount; levelIndex > 0;)
    {
        --levelIndex;
        if (searchShouldStop()) return;
        dagDuplicateClassLevel &stmap = frame.duplicateLevels[levelIndex];
        frame.aggregateMasks.reset(input.fragments.size() + 2);
        span<EdgeMaskAccumulator> stmapMaskList =
            frame.aggregateMasks.span().first(input.fragments.size());
        EdgeMaskAccumulator &maskC =
            frame.aggregateMasks[input.fragments.size()];
        EdgeMaskAccumulator &maskM =
            frame.aggregateMasks[input.fragments.size() + 1];
        vi &unrestrictedParentTotals = frame.unrestrictedParentTotals;
        vi &pairGenericBoundCache = frame.pairGenericBoundCache;
        unrestrictedParentTotals.clear();
        pairGenericBoundCache.clear();
        for (auto &entry : stmap.classes)
        {
            if (searchShouldStop()) return;
            dagDuplicateSet &ss = entry.duplicates;
            if (!ss.dead)
            {
            maskC.clear();
            const duplicateFragmentMaskList duplicateMasks =
                stmap.fragmentMasks(entry);
            
            for (const duplicateFragmentMaskEntry duplicateMask : duplicateMasks)
            {
                if (searchShouldStop()) return;
                maskC |= duplicateMask.mask;
                stmapMaskList[duplicateMask.fragment] |= duplicateMask.mask;
            }
            maskM |= maskC;
            int dupBondsMaxFrag = input.maxDupBonds(ss.size, stmapMaskList);

            int temp = fragSizeListMax[ss.size - 2] - 1;

            if (ss.size > 2) temp = max(temp, fragSizeListMax[ss.size - 3]);

            /// Initial branch-and-bound before the fragmentation step
            int earlySDP = max(dupBondsMaxFrag, temp);

            int earlyAIBound = static_cast<int>(totalBonds) -
                input.sumDupBonds - 1 - earlySDP;
            if (earlyAIBound < AI)
            {
                int matchingClassBound = numeric_limits<int>::min();
                if (usePairBound && ss.size == 2 && ss.list.size() >= 48)
                {
                    const int pairBoundLimit = static_cast<int>(totalBonds) -
                        input.sumDupBonds - 1 - AI;
                    if (dupBondsMaxFrag - 1 <= pairBoundLimit)
                    {
                        if (
                            dupBondsMaxFrag <= pairBoundLimit ||
                            (
                                matchingClassBound = input.maxDupBonds(
                                    ss.size,
                                    duplicateMasks
                                )
                            ) <= pairBoundLimit
                        )
                        {
                            continue;
                        }
                    }
                }
                auto pairBoundFiltersMatching = [&](validMatchings &matching)
                {
                    if (usePairBound)
                    {
                        const int pairBoundLimit = static_cast<int>(totalBonds) -
                            input.sumDupBonds - 1 - AI;
                        if (dupBondsMaxFrag - 1 <= pairBoundLimit)
                        {
                            bool targetedBoundsFit =
                                dupBondsMaxFrag <= pairBoundLimit;
                            if (!targetedBoundsFit)
                            {
                                if (
                                    matchingClassBound ==
                                    numeric_limits<int>::min()
                                )
                                {
                                    matchingClassBound = input.maxDupBonds(
                                        ss.size,
                                        duplicateMasks
                                    );
                                }
                                targetedBoundsFit =
                                    matchingClassBound <= pairBoundLimit;
                            }

                            if (targetedBoundsFit)
                            {
                                bool genericBoundFits =
                                    matching.maxFragSize == 2;
                                if (!genericBoundFits)
                                {
                                    if (unrestrictedParentTotals.empty())
                                    {
                                        buildUnrestrictedDupBondTotals(
                                            input,
                                            matching.maxFragSize - 1,
                                            unrestrictedParentTotals
                                        );
                                        pairGenericBoundCache.assign(
                                            input.fragments.size() *
                                                input.fragments.size(),
                                            numeric_limits<int>::min()
                                        );
                                    }
                                    const size_t firstFragment = min(
                                        matching.frag1,
                                        matching.frag2
                                    );
                                    const size_t secondFragment = max(
                                        matching.frag1,
                                        matching.frag2
                                    );
                                    int &genericRouteBound = pairGenericBoundCache[
                                        firstFragment * input.fragments.size() +
                                        secondFragment
                                    ];
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                                    if (searchTelemetryEnabled) [[unlikely]]
                                    {
                                        ++searchTelemetry.counters
                                            .pairBoundCacheLookups;
                                        if (
                                            genericRouteBound ==
                                            numeric_limits<int>::min()
                                        )
                                        {
                                            ++searchTelemetry.counters
                                                .pairBoundCacheMisses;
                                        }
                                        else
                                        {
                                            ++searchTelemetry.counters
                                                .pairBoundCacheHits;
                                        }
                                    }
#endif
                                    if (
                                        genericRouteBound ==
                                        numeric_limits<int>::min()
                                    )
                                    {
                                        genericRouteBound =
                                            matching.maxFragSize - 1 +
                                            pairSpecificGenericBound(
                                                input,
                                                matching,
                                                unrestrictedParentTotals
                                            );
                                    }
                                    genericBoundFits =
                                        genericRouteBound <= pairBoundLimit;
                                }
                                if (genericBoundFits) return true;
                            }
                        }
                    }
                    return false;
                };
                auto matchingVisitor = [&](validMatchings &matching)
                {
                    if constexpr (
                        equivalenceMode != matchingEquivalenceMode::none
                    )
                    {
                        if (pairBoundFiltersMatching(matching)) return true;
                    }
                    candidate.clearFragments();
                    fragmentAssemblyStateWithoutCanonisationWithWorkspace(
                        input,
                        matching,
                        entry.canonicalId,
                        candidate,
                        fragmentationWorkspace
                    );
                    if (searchShouldStop()) return false;

                    int sumDupBonds =
                        input.sumDupBonds + matching.maxFragSize - 1;
                    candidate.sumDupBonds = sumDupBonds;
                    int fragmentationCutoff = postFragmentationCutoff(
                        candidate,
                        maskC,
                        maskM,
                        fragmentationWorkspace.boundTotals
                    );
                    if (searchShouldStop()) return false;
                    const int candidateAIBound =
                        static_cast<int>(totalBonds) - sumDupBonds - 1 -
                        fragmentationCutoff;
                    if (candidateAIBound < AI)
                    {
#if defined(ASSEMBLYCPP_USE_OPENMP) || defined(ASSEMBLYCPP_USE_MPI)
                        if constexpr (allowDepthTwoDonation && !trackPath)
                        {
                            if (
                                parallelTaskScheduler != nullptr &&
                                parallelTaskScheduler
                                    ->depthTwoRefillRequested() &&
                                parallelTaskScheduler->tryEnqueueDepthTwo(
                                    candidate,
                                    candidateAIBound
                                )
                            ) [[unlikely]]
                            {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                                ++searchDepthTwoTasksSpawned;
#endif
                                return true;
                            }
                        }
#endif
                        vi &candidateKey = searchStorage.candidateKey;
                        if (!canoniseAssemblyStateAndBuildKey(
                            candidate,
                            candidateKey,
                            fragmentationWorkspace
                        )) return false;
                        // The default false specialization is deliberate:
                        // only this root state's immediate children may be
                        // donated, never descendants of an inline child.
                        if (!continueAssemblySearchWithWorkspace<
                            equivalenceMode,
                            trackPath,
                            false,
                            useSharedStates,
                            enableFragmentPairBlocks
                                ? fragmentPairBlockMode::enabled
                                : fragmentPairBlockMode::disabled
                        >(
                            dag,
                            candidate,
                            matching,
                            candidateKey,
                            sumDupBonds,
                            AI,
                            fragmentationWorkspace,
                            searchStorage
                        )) return false;
                    }
                    return true;
                };
                bool completed;
                if constexpr (
                    equivalenceMode ==
                    matchingEquivalenceMode::homogeneousPath
                )
                {
                    homogeneousPathEquivalentMatchings equivalentMatchings(
                        input,
                        ss,
                        fragmentationWorkspace.homogeneousPathEdgePositions
                    );
                    if (equivalentMatchings.enabled)
                    {
                        completed = ss.visitMatchingsInReverse(
                            [&](validMatchings &matching,
                                size_t firstOccurrence,
                                size_t secondOccurrence)
                            {
                                return equivalentMatchings.skip(
                                    matching,
                                    firstOccurrence,
                                    secondOccurrence
                                );
                            },
                            matchingVisitor
                        );
                    }
                    else
                    {
                        completed = ss.visitMatchingsInReverse(
                            matchingVisitor
                        );
                    }
                }
                else
                {
                    if constexpr (enableFragmentPairBlocks)
                    {
                        constexpr size_t minimumGroupedOccurrences = 8;
                        const int pairBoundLimit =
                            static_cast<int>(totalBonds) -
                            input.sumDupBonds - 1 - AI;
                        // Pathway output retains the legacy reverse-occurrence
                        // order. Larger one-word count-only classes can visit
                        // their existing fragment-contiguous runs as Cartesian
                        // blocks.
                        if (
                            usePairBound &&
                            ss.size > 2 &&
                            ss.list.size() >= minimumGroupedOccurrences &&
                            shouldVisitFragmentPairBlocks(
                                input,
                                ss,
                                duplicateMasks,
                                dupBondsMaxFrag,
                                pairBoundLimit,
                                matchingClassBound
                            )
                        )
                        {
                            // This callback is deliberately separate from the
                            // occurrence-based legacy filter above. Keeping
                            // the fragment-only path inside grouped traversal
                            // avoids adding an adapter call to every legacy
                            // occurrence pair.
                            auto pairBoundFiltersFragmentPair = [&]
                            (
                                int firstFragment,
                                int secondFragment,
                                int selectedSize
                            )
                            {
                                if (usePairBound)
                                {
                                    const int pairBoundLimit =
                                        static_cast<int>(totalBonds) -
                                        input.sumDupBonds - 1 - AI;
                                    if (
                                        dupBondsMaxFrag - 1 <= pairBoundLimit
                                    )
                                    {
                                        bool targetedBoundsFit =
                                            dupBondsMaxFrag <= pairBoundLimit;
                                        if (!targetedBoundsFit)
                                        {
                                            if (
                                                matchingClassBound ==
                                                numeric_limits<int>::min()
                                            )
                                            {
                                                matchingClassBound =
                                                    input.maxDupBonds(
                                                        ss.size,
                                                        duplicateMasks
                                                    );
                                            }
                                            targetedBoundsFit =
                                                matchingClassBound <=
                                                pairBoundLimit;
                                        }

                                        if (targetedBoundsFit)
                                        {
                                            bool genericBoundFits =
                                                selectedSize == 2;
                                            if (!genericBoundFits)
                                            {
                                                if (
                                                    unrestrictedParentTotals
                                                        .empty()
                                                )
                                                {
                                                    buildUnrestrictedDupBondTotals(
                                                        input,
                                                        selectedSize - 1,
                                                        unrestrictedParentTotals
                                                    );
                                                    pairGenericBoundCache.assign(
                                                        input.fragments.size() *
                                                            input.fragments.size(),
                                                        numeric_limits<int>::min()
                                                    );
                                                }
                                                const size_t cacheFirstFragment =
                                                    min(
                                                        firstFragment,
                                                        secondFragment
                                                    );
                                                const size_t cacheSecondFragment =
                                                    max(
                                                        firstFragment,
                                                        secondFragment
                                                    );
                                                int &genericRouteBound =
                                                    pairGenericBoundCache[
                                                        cacheFirstFragment *
                                                            input.fragments.size() +
                                                        cacheSecondFragment
                                                    ];
#ifdef ASSEMBLY_ENABLE_TELEMETRY
                                                if (
                                                    searchTelemetryEnabled
                                                ) [[unlikely]]
                                                {
                                                    ++searchTelemetry.counters
                                                        .pairBoundCacheLookups;
                                                    if (
                                                        genericRouteBound ==
                                                        numeric_limits<int>::min()
                                                    )
                                                    {
                                                        ++searchTelemetry.counters
                                                            .pairBoundCacheMisses;
                                                    }
                                                    else
                                                    {
                                                        ++searchTelemetry.counters
                                                            .pairBoundCacheHits;
                                                    }
                                                }
#endif
                                                if (
                                                    genericRouteBound ==
                                                    numeric_limits<int>::min()
                                                )
                                                {
                                                    genericRouteBound =
                                                        selectedSize - 1 +
                                                        pairSpecificGenericBound(
                                                            input,
                                                            selectedSize,
                                                            firstFragment,
                                                            secondFragment,
                                                            unrestrictedParentTotals
                                                        );
                                                }
                                                genericBoundFits =
                                                    genericRouteBound <=
                                                    pairBoundLimit;
                                            }
                                            if (genericBoundFits) return true;
                                        }
                                    }
                                }
                                return false;
                            };
                            completed = ss.visitMatchingsByFragmentPairInReverse(
                                pairBoundFiltersFragmentPair,
                                matchingVisitor
                            );
                        }
                        else
                        {
                            completed = ss.visitMatchingsInReverse(
                                [&](validMatchings &matching, size_t, size_t)
                                {
                                    return pairBoundFiltersMatching(matching);
                                },
                                matchingVisitor
                            );
                        }
                    }
                    else
                    {
                        completed = ss.visitMatchingsInReverse(
                            [&](validMatchings &matching, size_t, size_t)
                            {
                                return pairBoundFiltersMatching(matching);
                            },
                            matchingVisitor
                        );
                    }
                }
                if (!completed) return;
            }
            }
        }
    }
}

/**
 * @brief The recursive function that enumerates duplicates and generates assembly states on the first pass
 * of the assembly algorithm
 * 
 * @param input The input assembly state
 * @param AI The global minimum assembly index found
 * @param fragmentationWorkspace Buffers reused across the search
 */
template<matchingEquivalenceMode equivalenceMode, bool trackPath>
void initialRecursiveAssemblyWithWorkspaceImpl(
    vector<dagLevel> &dag,
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    recordImprovedAssemblyIndex<trackPath>(input, AI, searchStorage);
    if (searchShouldStop()) return;

    vector<initialDuplicateClassLevel> stmapVector;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    setSearchTelemetryPhase(SearchTelemetryPhase::initialEnumeration);
#endif
    const bool hasInitialMatchings = initialRecursiveEnumeration(
        input,
        stmapVector,
        searchStorage.duplicateClassIndex,
        dag
    );
    if (enumerationLimitReached || searchShouldStop()) return;

    if (!hasInitialMatchings) return;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    setSearchTelemetryPhase(SearchTelemetryPhase::assemblySearch);
#endif
    assemblyState candidate;
    candidate.reserveFragments(input.fragments.size() + 2);
    size_t branchLeaseBegin = 0;
    size_t branchLeaseEnd = 0;
    bool branchLeaseUsed = false;
    auto claimNextBranchLease = [&]() -> bool
    {
        size_t begin = sharedBranchLeaseCursor->load(std::memory_order_relaxed);
        while (true)
        {
            if (begin == std::numeric_limits<size_t>::max())
            {
                branchLeaseBegin = begin;
                branchLeaseEnd = begin;
                branchLeaseUsed = false;
                return false;
            }
            const size_t remaining =
                std::numeric_limits<size_t>::max() - begin;
            const size_t end = searchBranchLeaseSize > remaining
                ? std::numeric_limits<size_t>::max()
                : begin + searchBranchLeaseSize;
            if (sharedBranchLeaseCursor->compare_exchange_weak(
                begin,
                end,
                std::memory_order_relaxed,
                std::memory_order_relaxed
            ))
            {
                branchLeaseBegin = begin;
                branchLeaseEnd = end;
                branchLeaseUsed = false;
                return true;
            }
        }
    };
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        if (searchShouldStop()) return;
        initialDuplicateClassLevel &stmap = stmapVector[j];
        for (auto &entry : stmap.classes)
        {
            if (searchShouldStop()) return;
            initialDuplicateSet &ss = entry.duplicates;
            auto matchingVisitor = [&](validMatchings &matching)
            {
                const size_t branchOrdinal = searchRootBranchOrdinal++;
                if (sharedBranchLeaseCursor != nullptr)
                {
                    if (
                        branchOrdinal % searchRankPartitionCount !=
                        searchRankPartitionIndex
                    ) return true;
                    const size_t partitionOrdinal =
                        branchOrdinal / searchRankPartitionCount;
                    while (partitionOrdinal >= branchLeaseEnd)
                    {
                        if (!claimNextBranchLease()) return true;
                        if (partitionOrdinal < branchLeaseBegin) return true;
                    }
                    if (partitionOrdinal < branchLeaseBegin) return true;
                    if (!branchLeaseUsed)
                    {
                        ++searchBranchLeaseCount;
                        branchLeaseUsed = true;
                    }
                    ++searchBranchAssignmentCount;
                }
                if (sharedAssemblyIndex != nullptr)
                    AI = min(
                        AI,
                        sharedAssemblyIndex->load(std::memory_order_relaxed)
                    );
                candidate.clearFragments();
                fragmentAssemblyStateWithoutCanonisationWithWorkspace(
                    input,
                    matching,
                    entry.canonicalId,
                    candidate,
                    fragmentationWorkspace
                );
                if (searchShouldStop()) return false;
                int sumDupBonds =
                    input.sumDupBonds + matching.maxFragSize - 1;
                candidate.sumDupBonds = sumDupBonds;
                if (candidate.lowBoundAI() < AI)
                {
                    vi &candidateKey = searchStorage.candidateKey;
                    if (!canoniseAssemblyStateAndBuildKey(
                        candidate,
                        candidateKey,
                        fragmentationWorkspace
                    )) return false;
                    if (!continueAssemblySearchWithWorkspace<
                        equivalenceMode,
                        trackPath
                    >(
                        dag,
                        candidate,
                        matching,
                        candidateKey,
                        sumDupBonds,
                        AI,
                        fragmentationWorkspace,
                        searchStorage
                    )) return false;
                }
                return true;
            };
            bool completed;
            if constexpr (
                equivalenceMode ==
                matchingEquivalenceMode::homogeneousPath
            )
            {
                homogeneousPathEquivalentMatchings equivalentMatchings(
                    input,
                    ss,
                    fragmentationWorkspace.homogeneousPathEdgePositions
                );
                if (equivalentMatchings.enabled)
                {
                    completed = ss.visitMatchingsInReverse(
                        [&](validMatchings &matching,
                            size_t firstOccurrence,
                            size_t secondOccurrence)
                        {
                            return equivalentMatchings.skip(
                                matching,
                                firstOccurrence,
                                secondOccurrence
                            );
                        },
                        matchingVisitor
                    );
                }
                else
                {
                    completed = ss.visitMatchingsInReverse(matchingVisitor);
                }
            }
            else
            {
                completed = ss.visitMatchingsInReverse(matchingVisitor);
            }
            if (!completed) return;
        }
    }
}

template<bool trackPath>
void initialRecursiveAssemblyWithWorkspace(
    vector<dagLevel> &dag,
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    if (!fragmentationWorkspace.homogeneousPathEdgePositions.empty())
    {
        initialRecursiveAssemblyWithWorkspaceImpl<
            matchingEquivalenceMode::homogeneousPath,
            trackPath
        >(dag, input, AI, fragmentationWorkspace, searchStorage);
    }
    else
    {
        initialRecursiveAssemblyWithWorkspaceImpl<
            matchingEquivalenceMode::none,
            trackPath
        >(dag, input, AI, fragmentationWorkspace, searchStorage);
    }
}

template<bool trackPath>
bool runImprovedAssemblySearch(
    vector<dagLevel> &dag,
    assemblyState &root,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    vector<edgeL> &removedEdges,
    ofstream &ofs,
    assemblySearchStorage &searchStorage
)
{
    root.assemblyHashCalculator(searchStorage.candidateKey);
    static_cast<void>(searchStorage.states.consider(
        searchStorage.candidateKey,
        0
    ));

    initialRecursiveAssemblyWithWorkspace<trackPath>(
        dag,
        root,
        AI,
        fragmentationWorkspace,
        searchStorage
    );

#ifdef ASSEMBLY_ENABLE_TELEMETRY
    setSearchTelemetryPhase(SearchTelemetryPhase::output);
#endif
    // Append the numeric result even when the search stops at a limit.
    lastCalculatedAssemblyIndex = compensateDisjointAssemblyIndex(AI);
    ofs << lastCalculatedAssemblyIndex << '\n';
    if (runtimeLimitReached)
    {
        if (!suppressSearchOutput)
        {
#ifdef ASSEMBLYCPP_LIBRARY_BUILD
            if (verbose)
#endif
            cout << "status: runtime limit reached\n";
        }
        ofs << "status: runtime limit reached\n";
    }
    if (enumerationLimitReached)
    {
        if (!suppressSearchOutput)
        {
#ifdef ASSEMBLYCPP_LIBRARY_BUILD
            if (verbose)
#endif
            cout << "status: enumeration limit reached\n";
        }
        ofs << "status: enumeration limit reached\n";
    }

    if constexpr (trackPath)
    {
        return recoverPathway2(searchStorage.pathway->best, removedEdges);
    }
    return true;
}

/** Serialize the root occurrences once and retain only primitive job indices. */
bool buildRootJobDescriptors(
    assemblyState &root,
    vector<initialDuplicateClassLevel> &levels,
    SearchContext &context
)
{
    const size_t wordCount = EdgeMask::activeWordCount();
    size_t totalOccurrences = 0;
    for (const initialDuplicateClassLevel &level : levels)
    {
        for (const auto &entry : level.classes)
        {
            const size_t count = entry.duplicates.list.size();
            if (count > numeric_limits<size_t>::max() - totalOccurrences)
                throw length_error("root occurrences exceed capacity");
            totalOccurrences += count;
        }
    }
    if (
        wordCount != 0 &&
        totalOccurrences >
            numeric_limits<size_t>::max() / wordCount
    )
    {
        throw length_error("serialized root masks exceed capacity");
    }
    context.rootOccurrences.reserve(totalOccurrences);
    context.occurrenceWords.reserve(totalOccurrences * wordCount);

    for (int levelIndex = static_cast<int>(levels.size()) - 1;
         levelIndex >= 0;
         --levelIndex)
    {
        if (searchShouldStop()) return false;
        initialDuplicateClassLevel &level = levels[levelIndex];
        for (auto &entry : level.classes)
        {
            if (searchShouldStop()) return false;
            initialDuplicateSet &duplicates = entry.duplicates;
            if (
                duplicates.size >
                static_cast<size_t>(numeric_limits<uint32_t>::max()) ||
                duplicates.size >
                static_cast<size_t>(numeric_limits<int>::max())
            )
            {
                throw length_error("root duplicate size exceeds job capacity");
            }

            const size_t occurrenceBase = context.rootOccurrences.size();
            if (
                duplicates.list.size() >
                numeric_limits<size_t>::max() - occurrenceBase
            )
            {
                throw length_error("root occurrences exceed capacity");
            }
            if (
                wordCount != 0 &&
                duplicates.list.size() >
                    (numeric_limits<size_t>::max() -
                        context.occurrenceWords.size()) / wordCount
            )
            {
                throw length_error("serialized root masks exceed capacity");
            }
            for (const initialPotentialDuplicate &occurrence : duplicates.list)
            {
                context.rootOccurrences.push_back({
                    context.occurrenceWords.size(),
                    static_cast<int32_t>(occurrence.fragment)
                });
                for (size_t word = 0; word < wordCount; ++word)
                {
                    context.occurrenceWords.push_back(
                        occurrence.mask.activeWord(word)
                    );
                }
            }

            size_t firstOccurrence = 0;
            size_t secondOccurrence = 0;
            auto recordPair = [&](validMatchings &)
            {
                context.rootJobs.push_back({
                    occurrenceBase + firstOccurrence,
                    occurrenceBase + secondOccurrence,
                    static_cast<int32_t>(entry.canonicalId),
                    static_cast<uint32_t>(duplicates.size)
                });
                return true;
            };

            bool completed;
            if (!context.homogeneousPathEdgePositions.empty())
            {
                homogeneousPathEquivalentMatchings equivalentMatchings(
                    root,
                    duplicates,
                    context.homogeneousPathEdgePositions
                );
                completed = duplicates.visitMatchingsInReverse(
                    [&](validMatchings &matching,
                        size_t first,
                        size_t second)
                    {
                        firstOccurrence = first;
                        secondOccurrence = second;
                        return equivalentMatchings.enabled &&
                            equivalentMatchings.skip(matching, first, second);
                    },
                    recordPair
                );
            }
            else
            {
                completed = duplicates.visitMatchingsInReverse(
                    [&](validMatchings &, size_t first, size_t second)
                    {
                        firstOccurrence = first;
                        secondOccurrence = second;
                        return false;
                    },
                    recordPair
                );
            }
            if (!completed) return false;
        }
    }

    // Levels were visited largest-first, so the leading job is already the
    // strongest cheap incumbent candidate without a separate sorting pass.
    return true;
}

/**
 * Build the mask-free, read-only state shared by one process's workers.
 * Root-owned masks are explicitly released before the producer thread can
 * enter the worker pool and reconfigure its thread-local mask arena.
 */
void prepareParallelSearchContext(
    molGraph &mg,
    SearchContext &context,
    std::size_t localWorkerCount
)
{
    context.startedAt = clock();
    startTime = context.startedAt;
    searchStopPollCountdown = 0;
    searchStopInnerPollCountdown = 0;
    runtimeLimitReached = false;
    enumerationLimitReached = false;
    sharedCanonicalRegistry = nullptr;
    sharedAssemblyStates = nullptr;
    sharedAssemblyWorkerIndex = 0;
    context.canonicalRegistry.reset();
    context.sharedStates.reset();
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    intermediateMAs.clear();

    totalBonds = mg.totalBonds;
    disjointFragments = mg.disjointFragments();
    vector<edgeL> removedEdges;
    targetMolecule = preprocessWriteback(mg, removedEdges);
    univEdgeList = targetMolecule.writeEdgeList();
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);

    // Release the persistent mask before configure() clears this TLS arena.
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(univEdgeList.size());
    std::construct_at(std::addressof(allEdges));
    AtomMask::configure(targetMolecule.mg.size());

    context.homogeneousPathEdgePositions.clear();
    configureHomogeneousPathEdgePositions(
        context.homogeneousPathEdgePositions
    );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    configureSearchTelemetryGraph(
        targetMolecule.mg.size(),
        univEdgeList.size(),
        EdgeMask::activeWordCount(),
        ufdsMaskWorkspace::decompositionCacheEligible(univEdgeList.size())
    );
    setSearchTelemetryPhase(SearchTelemetryPhase::initialEnumeration);
#endif

    allEdges.set();
    assemblyState root;
    root.appendFragment(
        allEdges,
        static_cast<int>(univEdgeList.size()),
        unknownCanonicalId,
        false
    );
    context.rootAssemblyIndex = root.AI();

    vector<initialDuplicateClassLevel> levels;
    duplicateClassIndexWorkspace classIndex;
    const bool hasInitialMatchings = initialRecursiveEnumeration(
        root,
        levels,
        classIndex,
        context.dag
    );
    if (
        hasInitialMatchings &&
        !enumerationLimitReached &&
        !searchShouldStop()
    )
    {
        static_cast<void>(buildRootJobDescriptors(root, levels, context));
    }
    // Locking cannot repay its fixed setup cost on short or low-reuse graphs.
    // Small homogeneous paths also finish almost entirely in their dedicated
    // root-equivalence pass; larger frontiers leave enough work for L2 reuse.
    const bool useSharedReuse =
        localWorkerCount > 1 &&
        univEdgeList.size() > 30 &&
        graphHashMap.size() <= bitsetHashTable.size() / 2 &&
        (
            context.homogeneousPathEdgePositions.empty() ||
            context.rootJobs.size() >= 10000
        );
    if (useSharedReuse)
    {
        context.canonicalRegistry =
            make_unique<sharedCanonicalIdRegistry>(graphHashMap.size());
        context.sharedStates =
            make_unique<sharedAssemblyTranspositionTable>(localWorkerCount);
    }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    setSearchTelemetryPhase(SearchTelemetryPhase::assemblySearch);
#endif

    context.bondCount = totalBonds;
    context.componentCount = disjointFragments;
    context.enumerationLimit = enumerationLimitReached;

    // These caches contain no active-word masks and can seed worker-local
    // canonical IDs. Move them only after DAG conversion and job serialization.
    context.canonicalSeed.graphHashes = std::move(graphHashMap);
    context.canonicalSeed.atomInterner = std::move(treeCanonAtomInterner);
    context.canonicalSeed.leafInterner = std::move(treeCanonLeafInterner);
    context.canonicalSeed.treeInterner = std::move(treeCanonInterner);

    // Destroy every producer-owned wide mask before this same thread becomes
    // OpenMP worker zero and calls EdgeMask::configure again.
    levels.clear();
    root.clearFragments();
    bitsetHashTable.clear();
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(univEdgeList.size());
    std::construct_at(std::addressof(allEdges));

    context.processedMolecule = std::move(targetMolecule);
    context.universeEdges = std::move(univEdgeList);
}

/** Configure only thread-local state from an immutable, mask-free context. */
void configureParallelWorker(
    const SearchContext &context,
    std::size_t workerIndex
)
{
    if (context.sharedStates != nullptr)
    {
        sharedCanonicalRegistry = context.canonicalRegistry.get();
        sharedAssemblyStates = context.sharedStates.get();
        sharedAssemblyWorkerIndex = workerIndex;
    }
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(context.universeEdges.size());
    std::construct_at(std::addressof(allEdges));
    AtomMask::configure(context.processedMolecule.mg.size());

    startTime = context.startedAt;
    searchStopPollCountdown = 0;
    searchStopInnerPollCountdown = 0;
    runtimeLimitReached = false;
    enumerationLimitReached = context.enumerationLimit;
    totalBonds = context.bondCount;
    disjointFragments = context.componentCount;
    targetMolecule = context.processedMolecule;
    univEdgeList = context.universeEdges;

    treeCanonAtomInterner = context.canonicalSeed.atomInterner;
    treeCanonLeafInterner = context.canonicalSeed.leafInterner;
    treeCanonInterner = context.canonicalSeed.treeInterner;
    graphHashMap = context.canonicalSeed.graphHashes;
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    configureSearchTelemetryGraph(
        targetMolecule.mg.size(),
        univEdgeList.size(),
        EdgeMask::activeWordCount(),
        ufdsMaskWorkspace::decompositionCacheEligible(univEdgeList.size())
    );
#endif
}

/** Release all worker-owned masks while still on their allocating thread. */
void clearParallelWorkerMasks()
{
    bitsetHashTable.clear();
    allEdges.reset();
    if (sharedAssemblyStates != nullptr)
    {
        sharedCanonicalRegistry = nullptr;
        sharedAssemblyStates = nullptr;
        sharedAssemblyWorkerIndex = 0;
    }
}

EdgeMask reconstructRootOccurrence(
    const SearchContext &context,
    size_t occurrenceIndex
)
{
    if (occurrenceIndex >= context.rootOccurrences.size())
        throw logic_error("root job occurrence is out of range");
    const rootOccurrenceDescriptor &occurrence =
        context.rootOccurrences[occurrenceIndex];
    const size_t wordCount = EdgeMask::activeWordCount();
    if (
        occurrence.wordOffset > context.occurrenceWords.size() ||
        wordCount > context.occurrenceWords.size() - occurrence.wordOffset
    )
    {
        throw logic_error("serialized root mask is incomplete");
    }
    return EdgeMask::fromActiveWords(
        context.occurrenceWords.data() + occurrence.wordOffset
    );
}

template<matchingEquivalenceMode equivalenceMode, bool useSharedStates>
bool runParallelRootJobImpl(
    const SearchContext &context,
    size_t jobIndex,
    WorkerContext &worker
)
{
    if (jobIndex >= context.rootJobs.size())
        throw logic_error("root job index is out of range");
    const rootJobDescriptor &job = context.rootJobs[jobIndex];
    const rootOccurrenceDescriptor &firstOccurrence =
        context.rootOccurrences[job.firstOccurrence];
    const rootOccurrenceDescriptor &secondOccurrence =
        context.rootOccurrences[job.secondOccurrence];
    EdgeMask first = reconstructRootOccurrence(context, job.firstOccurrence);
    EdgeMask second = reconstructRootOccurrence(context, job.secondOccurrence);
    validMatchings matching(
        first,
        second,
        firstOccurrence.fragment,
        secondOccurrence.fragment,
        static_cast<int>(job.duplicateSize)
    );

    if (sharedAssemblyIndex != nullptr)
    {
        worker.assemblyIndex = min(
            worker.assemblyIndex,
            sharedAssemblyIndex->load(std::memory_order_relaxed)
        );
    }
    worker.candidate.clearFragments();
    fragmentAssemblyStateWithoutCanonisationWithWorkspace(
        worker.root,
        matching,
        job.canonicalId,
        worker.candidate,
        worker.fragmentation
    );
    if (searchShouldStop()) return false;

    const int sumDupBonds = matching.maxFragSize - 1;
    worker.candidate.sumDupBonds = sumDupBonds;
    if (worker.candidate.lowBoundAI() >= worker.assemblyIndex) return true;

    vi &candidateKey = worker.search.candidateKey;
    if (!canoniseAssemblyStateAndBuildKey(
        worker.candidate,
        candidateKey,
        worker.fragmentation
    )) return false;
    // The worker L1 removes repeated local root children. The exact shared L2
    // linearizes equivalent canonical children reached by different workers,
    // so only their first/best encounter proceeds into recursive search.
    return continueAssemblySearchWithWorkspace<
        equivalenceMode,
        false,
        true,
        useSharedStates
    >(
        context.dag,
        worker.candidate,
        matching,
        candidateKey,
        sumDupBonds,
        worker.assemblyIndex,
        worker.fragmentation,
        worker.search
    );
}

template<bool useSharedStates>
bool runParallelRootJob(
    const SearchContext &context,
    size_t jobIndex,
    WorkerContext &worker
)
{
    if (!worker.fragmentation.homogeneousPathEdgePositions.empty())
    {
        return runParallelRootJobImpl<
            matchingEquivalenceMode::homogeneousPath,
            useSharedStates
        >(context, jobIndex, worker);
    }
    return runParallelRootJobImpl<
        matchingEquivalenceMode::none,
        useSharedStates
    >(context, jobIndex, worker);
}

/** Evaluate the strongest one-step branch before workers enter the queue. */
void warmStartParallelIncumbent(
    const SearchContext &context,
    size_t jobIndex,
    WorkerContext &worker
)
{
    if (jobIndex >= context.rootJobs.size()) return;
    const rootJobDescriptor &job = context.rootJobs[jobIndex];
    const rootOccurrenceDescriptor &firstOccurrence =
        context.rootOccurrences[job.firstOccurrence];
    const rootOccurrenceDescriptor &secondOccurrence =
        context.rootOccurrences[job.secondOccurrence];
    EdgeMask first = reconstructRootOccurrence(context, job.firstOccurrence);
    EdgeMask second = reconstructRootOccurrence(context, job.secondOccurrence);
    validMatchings matching(
        first,
        second,
        firstOccurrence.fragment,
        secondOccurrence.fragment,
        static_cast<int>(job.duplicateSize)
    );
    worker.candidate.clearFragments();
    fragmentAssemblyStateWithoutCanonisationWithWorkspace(
        worker.root,
        matching,
        job.canonicalId,
        worker.candidate,
        worker.fragmentation
    );
    worker.candidate.sumDupBonds = matching.maxFragSize - 1;
    recordImprovedAssemblyIndex<false>(
        worker.candidate,
        worker.assemblyIndex,
        worker.search
    );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    ++searchWarmStartBranches;
#endif
}

void reconstructDepthTwoTask(
    const parallelDepthTwoTaskDescriptor &task,
    assemblyState &state
)
{
    state.clearFragments();
    state.reserveFragments(task.fragments.size());
    const size_t wordCount = EdgeMask::activeWordCount();
    for (const parallelTaskFragmentDescriptor &fragment : task.fragments)
    {
        if (
            fragment.wordOffset > task.fragmentWords.size() ||
            wordCount > task.fragmentWords.size() - fragment.wordOffset
        ) throw logic_error("serialized depth-two task mask is incomplete");
        EdgeMask mask = EdgeMask::fromActiveWords(
            task.fragmentWords.data() + fragment.wordOffset
        );
        state.appendFragment(
            mask,
            static_cast<int>(fragment.edgeCount),
            unknownCanonicalId,
            fragment.connected
        );
    }
    state.sumDupBonds = task.sumDupBonds;
}

template<matchingEquivalenceMode equivalenceMode, bool useSharedStates>
bool runParallelDepthTwoJobImpl(
    const SearchContext &context,
    const parallelDepthTwoTaskDescriptor &task,
    WorkerContext &worker
)
{
    if (sharedAssemblyIndex != nullptr)
    {
        worker.assemblyIndex = min(
            worker.assemblyIndex,
            sharedAssemblyIndex->load(std::memory_order_relaxed)
        );
    }
    if (task.lowerBoundAssemblyIndex >= worker.assemblyIndex) return true;
    // The previous root may have ended on a bound-pruned raw fragmentation.
    // Its deferred cache binding does not describe this transferred state.
    worker.fragmentation.beginFragmentation();
    reconstructDepthTwoTask(task, worker.candidate);
    if (searchShouldStop()) return false;

    vi &candidateKey = worker.search.candidateKey;
    if (!canoniseAssemblyStateAndBuildKey(
        worker.candidate,
        candidateKey,
        worker.fragmentation
    )) return false;
    return continueCanonicalAssemblySearchWithWorkspace<
        equivalenceMode,
        false,
        false,
        useSharedStates
    >(
        context.dag,
        worker.candidate,
        candidateKey,
        worker.candidate.sumDupBonds,
        worker.assemblyIndex,
        worker.fragmentation,
        worker.search
    );
}

template<bool useSharedStates>
bool runParallelDepthTwoJob(
    const SearchContext &context,
    const parallelDepthTwoTaskDescriptor &task,
    WorkerContext &worker
)
{
    if (!worker.fragmentation.homogeneousPathEdgePositions.empty())
    {
        return runParallelDepthTwoJobImpl<
            matchingEquivalenceMode::homogeneousPath,
            useSharedStates
        >(context, task, worker);
    }
    return runParallelDepthTwoJobImpl<
        matchingEquivalenceMode::none,
        useSharedStates
    >(context, task, worker);
}

/** Dynamically lease roots and consume adaptively exposed depth-two work. */
template<bool useSharedStates>
void runParallelRootJobs(
    const SearchContext &context,
    WorkerContext &worker,
    ParallelTaskScheduler &scheduler
)
{
    searchRootBranchOrdinal = context.rootJobs.size();
    while (!searchShouldStop())
    {
        size_t begin = 0;
        size_t end = 0;
        if (scheduler.claimRootLease(begin, end))
        {
            ++searchBranchLeaseCount;
            for (size_t partitionOrdinal = begin;
                 partitionOrdinal < end;
                 ++partitionOrdinal)
            {
                const size_t jobIndex = searchRankPartitionIndex +
                    partitionOrdinal * searchRankPartitionCount;
                ++searchBranchAssignmentCount;
                bool completed = false;
                try
                {
                    completed = runParallelRootJob<useSharedStates>(
                        context,
                        jobIndex,
                        worker
                    );
                }
                catch (...)
                {
                    scheduler.completeRootLease(partitionOrdinal - begin);
                    throw;
                }
                if (!completed)
                {
                    scheduler.completeRootLease(
                        partitionOrdinal - begin + 1
                    );
                    return;
                }
            }
            scheduler.completeRootLease(end - begin);
            continue;
        }

        if (!scheduler.depthTwoTasksEnabled()) return;

        // Root work remains the normal scheduling frontier. Depth-two tasks
        // are consumed only after no root lease is immediately claimable.
        parallelDepthTwoTaskDescriptor task;
        const ParallelTaskScheduler::WorkAvailability availability =
            scheduler.nextWork(task);
        if (
            availability ==
                ParallelTaskScheduler::WorkAvailability::depthTwo
        )
        {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            ++searchDepthTwoTasksExecuted;
#endif
            bool completed = false;
            try
            {
                completed = runParallelDepthTwoJob<useSharedStates>(
                    context,
                    task,
                    worker
                );
            }
            catch (...)
            {
                scheduler.completeDepthTwoJob();
                throw;
            }
            scheduler.completeDepthTwoJob();
            if (!completed) return;
            continue;
        }

        if (
            availability == ParallelTaskScheduler::WorkAvailability::complete
        ) return;
        scheduler.waitForWork();
    }
}

/**
 * @brief Function that calls the recursive assembly function
 * 
 * @param mg The target molGraph
 * @param ofs The output file
 * @return false only when a requested pathway could not be written.
 */
bool improvedBnB(molGraph &mg, ofstream &ofs)
{
    startTime = clock();
    searchStopPollCountdown = 0;
    searchStopInnerPollCountdown = 0;
    runtimeLimitReached = false;
    enumerationLimitReached = false;
    sharedCanonicalRegistry = nullptr;
    sharedAssemblyStates = nullptr;
    sharedAssemblyWorkerIndex = 0;
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    intermediateMAs.clear();
    searchRootBranchOrdinal = 0;
    searchBranchLeaseCount = 0;
    searchBranchAssignmentCount = 0;
    totalBonds = mg.totalBonds;
    originalEdgeList = mg.writeEdgeList();
    disjointFragments = mg.disjointFragments();
    originalMolecule = mg;
    vector<edgeL> removedEdges;
    targetMolecule = preprocessWriteback(mg, removedEdges);
    univEdgeList = targetMolecule.writeEdgeList();
    prepareCanonicalisationGraph(targetMolecule, univEdgeList);
    // End the persistent mask's lifetime under its old representation before
    // changing the domain width, then construct its new representation. Other
    // mask-owning globals were cleared above.
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(univEdgeList.size());
    std::construct_at(std::addressof(allEdges));
    AtomMask::configure(targetMolecule.mg.size());
    ufdsMaskWorkspace fragmentationWorkspace(
        targetMolecule.mg.size(),
        univEdgeList.size()
    );
    configureHomogeneousPathEdgePositions(
        fragmentationWorkspace.homogeneousPathEdgePositions
    );
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    configureSearchTelemetryGraph(
        targetMolecule.mg.size(),
        univEdgeList.size(),
        EdgeMask::activeWordCount(),
        fragmentationWorkspace.reuseResidualDecompositions
    );
#endif
    for (size_t i = 0; i < univEdgeList.size(); i++) allEdges.set(i);
    assemblyState as;
    as.appendFragment(
        allEdges,
        static_cast<int>(univEdgeList.size()),
        unknownCanonicalId,
        false
    );
    vector<dagLevel> dag;
    int AI = std::numeric_limits<int>::max();
    if (isPathway)
    {
        assemblyPathWitness pathwayWitness;
        assemblySearchStorage searchStorage(&pathwayWitness);
        return runImprovedAssemblySearch<true>(
            dag,
            as,
            AI,
            fragmentationWorkspace,
            removedEdges,
            ofs,
            searchStorage
        );
    }

    assemblySearchStorage searchStorage;
    return runImprovedAssemblySearch<false>(
        dag,
        as,
        AI,
        fragmentationWorkspace,
        removedEdges,
        ofs,
        searchStorage
    );
}
