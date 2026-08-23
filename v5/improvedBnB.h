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
    duplicateClassIndexWorkspace &classIndex
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
                    incidentEdges
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
    convertDag(tempDag);
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
    assemblyState &_target,
    vector<dagDuplicateClassLevel> &stmapVector,
    vector<vector<EdgeMask>> &targetMasks,
    duplicateClassIndexWorkspace &classIndex
)
{
    if (searchShouldStop()) return 0;
    int ordinal = std::numeric_limits<int>::max();
    if (_target.fragments.front().canonicalId >= 0)
        ordinal = _target.fragments.front().canonicalId;
    vector<assemblyFragment> &fragments = _target.fragments;
    size_t currSize = 1;
    
    stmapVector.resize(1);
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
                    m,
                    stmapVector[0],
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
    stmapVector[0].seal();
    bool active = 1, overweight = 0, last = 0;
    while (active)
    {
        if (searchShouldStop()) return 0;
        vector<EdgeMask> targetMask(fragments.size(), 0);
        active = 0;
        stmapVector.emplace_back();
        classIndex.beginLevel();
        dagDuplicateClassLevel &stmap =
            stmapVector[stmapVector.size() - 2];
        for (auto &entry : stmap.classes)
        {
            if (searchShouldStop()) return 0;
            dagDuplicateSet &ss = entry.duplicates;
            if (ss.isValid())
            {
                if (searchShouldStop()) return 0;
                active |= dagDuplicateGenerator(
                    ss,
                    stmapVector.back(),
                    classIndex,
                    targetMask,
                    fragments,
                    ordinal,
                    overweight,
                    last
                );
                if (searchShouldStop()) return 0;
            }
        }
        stmapVector.back().seal();
        if (overweight) last = 1;
        targetMasks.push_back(targetMask);
        currSize++;
    }
    if (stmapVector.back().size() == 0) stmapVector.pop_back();
    for (size_t i = 0; i < fragments.size(); i++)
    {
        if (searchShouldStop()) return 0;
        assemblyFragment &fragment = fragments[i];
        fragment.retainEdges(targetMasks[0][i]);
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
int postFragmentationCutoff(
    const assemblyState &target,
    const EdgeMask &matchMask,
    const EdgeMask &maxFragMask,
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
                    fragment.mask.intersectionCount(matchMask)
                );
            const int maximalEdges = i == 0
                ? maxFragSize
                : static_cast<int>(
                    fragment.mask.intersectionCount(maxFragMask)
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

struct assemblyPathWitness
{
    vector<assemblyPathStep> current;
    vector<assemblyPathStep> best;
};

struct assemblySearchStorage
{
    assemblyTranspositionTable states;
    assemblyPathWitness *pathway = nullptr;
    duplicateClassIndexWorkspace duplicateClassIndex;
    vi candidateKey;

    explicit assemblySearchStorage(assemblyPathWitness *_pathway = nullptr):
        states(1024),
        pathway(_pathway)
    {
        // One search-wide scratch key is safe because every table operation
        // finishes before the synchronous recursive call can reuse it.
        candidateKey.reserve(univEdgeList.size() + 1);
        if (pathway != nullptr)
        {
            pathway->current.reserve(univEdgeList.size());
            pathway->best.reserve(univEdgeList.size());
        }
    }
};

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

template<matchingEquivalenceMode equivalenceMode, bool trackPath>
void dagRecursiveAssemblyWithWorkspaceImpl(
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

template<matchingEquivalenceMode equivalenceMode, bool trackPath>
bool continueAssemblySearchWithWorkspace(
    assemblyState &candidate,
    validMatchings &matching,
    span<const int> candidateKey,
    int sumDupBonds,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    if (searchShouldStop()) return false;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
        ++searchTelemetry.counters.assemblyCacheLookups;
#endif

    const assemblyTranspositionTable::result result =
        searchStorage.states.consider(candidateKey, sumDupBonds);

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
        searchStorage.pathway->current.push_back(
            assemblyPathStep{matching.first, matching.second}
        );
    }
    dagRecursiveAssemblyWithWorkspaceImpl<equivalenceMode, trackPath>(
        candidate,
        AI,
        fragmentationWorkspace,
        searchStorage
    );
    if constexpr (trackPath) searchStorage.pathway->current.pop_back();
    return !searchShouldStop();
}

/**
 * @brief The recursive function that enumerates duplicates and generates assembly states on all but the first pass
 * of the assembly algorithm
 * 
 * @param input The input assembly state
 * @param AI The global minimum assembly index found
 * @param fragmentationWorkspace Buffers reused across the search
 */
template<matchingEquivalenceMode equivalenceMode, bool trackPath>
void dagRecursiveAssemblyWithWorkspaceImpl(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    // Count-only pair pruning pays for itself once residual splitting is no
    // longer dominated by the fixed cost of the bound lookup.
    constexpr size_t pairBoundMinimumMoleculeEdges = 27;
    const bool usePairBound =
        fragmentationWorkspace.edgeCount >= pairBoundMinimumMoleculeEdges;
    recordImprovedAssemblyIndex<trackPath>(input, AI, searchStorage);
    if (searchShouldStop()) return;

    vector<dagDuplicateClassLevel> stmapVector;
    vector<vector<EdgeMask> > targetMasks;
    int maxFragSize = dagRecursiveEnumeration(
        input,
        stmapVector,
        targetMasks,
        searchStorage.duplicateClassIndex
    );

    if (searchShouldStop()) return;

    /// Find the fragment-size-specific AI lower bounds
    vi fragSizeListMax;
    input.maxDupBondsPrefix(fragSizeListMax, maxFragSize, targetMasks);
    if (searchShouldStop()) return;

    /// Begin iterating through the enumerated duplicatable fragments
    assemblyState candidate;
    candidate.reserveFragments(input.fragments.size() + 2);
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        if (searchShouldStop()) return;
        dagDuplicateClassLevel &stmap = stmapVector[j];
        vector<EdgeMask> stmapMaskList(input.fragments.size(), 0);
        vi unrestrictedParentTotals;
        vi pairGenericBoundCache;
        EdgeMask maskM = 0;
        for (auto &entry : stmap.classes)
        {
            if (searchShouldStop()) return;
            dagDuplicateSet &ss = entry.duplicates;
            if (!ss.dead)
            {
            EdgeMask maskC = 0;
            
            for (size_t i = 0; i < ss.maskList.size(); i++)
            {
                if (searchShouldStop()) return;
                maskC |= ss.maskList[i];
                stmapMaskList[i] |= ss.maskList[i];
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
                                    ss.maskList
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
                                        ss.maskList
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
                    completed = ss.visitMatchingsInReverse(
                        [&](validMatchings &matching, size_t, size_t)
                        {
                            return pairBoundFiltersMatching(matching);
                        },
                        matchingVisitor
                    );
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
        searchStorage.duplicateClassIndex
    );
    if (enumerationLimitReached || searchShouldStop()) return;

    if (!hasInitialMatchings) return;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    setSearchTelemetryPhase(SearchTelemetryPhase::assemblySearch);
#endif
    assemblyState candidate;
    candidate.reserveFragments(input.fragments.size() + 2);
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
                const size_t branchOrdinal = searchShardBranchOrdinal++;
                if (
                    searchShardCount > 1 &&
                    branchOrdinal % searchShardCount != searchShardIndex
                ) return true;
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
        >(input, AI, fragmentationWorkspace, searchStorage);
    }
    else
    {
        initialRecursiveAssemblyWithWorkspaceImpl<
            matchingEquivalenceMode::none,
            trackPath
        >(input, AI, fragmentationWorkspace, searchStorage);
    }
}

template<bool trackPath>
bool runImprovedAssemblySearch(
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
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    DAG.clear();
    intermediateMAs.clear();
    searchShardBranchOrdinal = 0;
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
    int AI = std::numeric_limits<int>::max();
    if (isPathway)
    {
        assemblyPathWitness pathwayWitness;
        assemblySearchStorage searchStorage(&pathwayWitness);
        return runImprovedAssemblySearch<true>(
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
        as,
        AI,
        fragmentationWorkspace,
        removedEdges,
        ofs,
        searchStorage
    );
}
