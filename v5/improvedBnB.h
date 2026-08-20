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
    vector<EdgeMask> &masks = _target.masks;
    bool alive = 0;
    size_t currSize = 1;
    vector<initialPotentialDuplicate> prevML;
    vector<int> rootNodeIndices(univEdgeList.size(), -1);
    size_t retainedStateCount = 0;

    // Retain the one-edge DAG states first so they count toward ENUM_MAX too.
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (searchShouldStop()) return false;
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStopPeriodically()) return false;
            if (masks[i][j] != 0)
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
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (searchShouldStop()) return false;
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStopPeriodically()) return false;
            if (masks[i][j] == 0) continue;

            if (rootNodeIndices[j] < 0)
                throw logic_error("initial DAG root was not retained");
            initialPotentialDuplicate m(
                j,
                masks[i],
                i,
                rootNodeIndices[j]
            );
            if (!m.generateDAG(prevML, retainedStateCount, tempDag))
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
                masks.size()
            ).insert(m);
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
                    tempDag
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
    const auto ordinalEntry = bitsetHashTable.find(_target.masks.front());
    if (ordinalEntry != bitsetHashTable.end()) ordinal = ordinalEntry->second.first;
    vector<EdgeMask> &masks = _target.masks;
    size_t currSize = 1;
    
    stmapVector.resize(1);
    classIndex.beginLevel();
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (searchShouldStop()) return 0;
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStopPeriodically()) return 0;
            if (masks[i][j] != 0)
            {
                EdgeMask b = 0; b.set(j);
                potentialDuplicate m(std::move(b), i, j);
                dagGenerate(
                    m,
                    stmapVector[0],
                    classIndex,
                    masks[i],
                    currSize,
                    ordinal,
                    masks.size()
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
        vector<EdgeMask> targetMask(masks.size(), 0);
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
                    masks,
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
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (searchShouldStop()) return 0;
        masks[i] &= targetMasks[0][i];
        _target.edgeCounts[i] = static_cast<int>(masks[i].count());
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
    const int maxFragSize = target.edgeCounts[0];
    const bool useTargetedBounds = target.masks.size() >= 2;
    int matchDB = 0;
    int maxFragDB = 0;
    boundTotals.assign(max(maxFragSize - 2, 0), 0);

    for (size_t i = 0; i < target.masks.size(); i++)
    {
        const int edgeCount = target.edgeCounts[i];
        if (useTargetedBounds)
        {
            const int matchingEdges = i == 0
                ? maxFragSize
                : static_cast<int>(
                    target.masks[i].intersectionCount(matchMask)
                );
            const int maximalEdges = i == 0
                ? maxFragSize
                : static_cast<int>(
                    target.masks[i].intersectionCount(maxFragMask)
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
    for (const int edgeCount : target.edgeCounts)
    {
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
        const int parentEdges = target.edgeCounts[matching.frag1];
        total += (parentEdges - 2 * selectedSize) / 2 - parentEdges / 2;
    }
    else
    {
        const int firstParentEdges = target.edgeCounts[matching.frag1];
        const int secondParentEdges = target.edgeCounts[matching.frag2];
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
            const int parentEdges = target.edgeCounts[matching.frag1];
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
            const int firstParentEdges = target.edgeCounts[matching.frag1];
            const int secondParentEdges = target.edgeCounts[matching.frag2];
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

struct assemblyStateKeyHash
{
    size_t operator()(const vi &key) const
    {
        size_t seed = key.size();
        for (const int value : key)
        {
            seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct pathwayTransposition
{
    int bestSumDupBonds;
    assemblyPath *path;
};

constexpr std::pmr::pool_options assemblyStatePoolOptions()
{
    return {
        .max_blocks_per_chunk = 4096,
        .largest_required_pool_block = 128
    };
}

struct scoreOnlyAssemblySearchStorage
{
    // Deliberately keep the no-path table to the state key and duplicate score.
    std::pmr::unsynchronized_pool_resource statePool{
        assemblyStatePoolOptions()
    };
    std::pmr::unordered_map<vi, int, assemblyStateKeyHash> states{&statePool};
};

struct pathwayAssemblySearchStorage
{
    static_assert(std::is_trivially_destructible_v<assemblyPath>);

    std::pmr::unsynchronized_pool_resource statePool{
        assemblyStatePoolOptions()
    };
    // Parent-link records have stable addresses and are released as one arena.
    std::pmr::monotonic_buffer_resource pathArena;
    std::pmr::unordered_map<vi, pathwayTransposition, assemblyStateKeyHash>
        states{&statePool};

    assemblyPath *createPath(
        int retainedFragmentClass,
        assemblyPath *parent,
        unsigned short match,
        unsigned short duplicate
    )
    {
        void *memory = pathArena.allocate(
            sizeof(assemblyPath),
            alignof(assemblyPath)
        );
        return std::construct_at(
            static_cast<assemblyPath *>(memory),
            assemblyPath{
                retainedFragmentClass,
                match,
                duplicate,
                parent
            }
        );
    }
};

struct assemblySearchStorage
{
    scoreOnlyAssemblySearchStorage *scoreOnly = nullptr;
    pathwayAssemblySearchStorage *pathway = nullptr;
    duplicateClassIndexWorkspace duplicateClassIndex;

    bool tracksPath() const noexcept
    {
        return pathway != nullptr;
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

        fragmentStarts.resize(input.masks.size());
        for (size_t fragment = 0; fragment < input.masks.size(); fragment++)
        {
            if (!intervalStart(
                input.masks[fragment],
                input.edgeCounts[fragment],
                fragmentStarts[fragment]
            )) return;
        }

        occurrenceStarts.resize(duplicates.list.size());
        vector<vector<int>> startsByFragment(input.masks.size());
        for (size_t occurrence = 0;
             occurrence < duplicates.list.size();
             occurrence++)
        {
            const auto &candidate = duplicates.list[occurrence];
            if (
                candidate.fragment < 0 ||
                static_cast<size_t>(candidate.fragment) >= input.masks.size() ||
                !input.masks[candidate.fragment].contains(candidate.mask) ||
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

            const int parentEdges = input.edgeCounts[fragment];
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
                const int parentEdges = input.edgeCounts[fragment];
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

void dagRecursiveAssemblyWithWorkspace(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
);

template<matchingEquivalenceMode equivalenceMode>
void dagRecursiveAssemblyWithWorkspaceImpl(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
);

void recordImprovedAssemblyIndex(
    assemblyState &input,
    int &AI,
    const assemblySearchStorage &searchStorage
)
{
    const int candidate = input.AI();
    if (candidate >= AI) return;

    AI = candidate;
    if (searchStorage.tracksPath()) minAssemblyPath = input.apPtr;
    const unsigned long long time = elapsedClockTicks();
    cout << "time: " << time << " min AI found so far: " << AI << '\n';
    if (writeIntermediateMAs) intermediateMAs.emplace_back(time, AI);
}

void setAssemblyPathStep(
    assemblyPath &path,
    assemblyPath *parent,
    validMatchings &matching
)
{
    path.parent = parent;
    path.match = bitsetHashTable.find(matching.first)->second.second;
    path.duplicate = bitsetHashTable.find(matching.second)->second.second;
}

template<matchingEquivalenceMode equivalenceMode>
bool continueAssemblySearchWithWorkspace(
    assemblyState &parent,
    assemblyState &candidate,
    validMatchings &matching,
    vi candidateKey,
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

    bool inserted = false;
    int *bestSumDupBonds = nullptr;
    assemblyPath *existingPath = nullptr;

    if (!searchStorage.tracksPath())
    {
        auto result = searchStorage.scoreOnly->states.try_emplace(
            std::move(candidateKey),
            sumDupBonds
        );
        inserted = result.second;
        bestSumDupBonds = &result.first->second;
    }
    else
    {
        pathwayAssemblySearchStorage &pathwayStorage = *searchStorage.pathway;
        const int retainedFragmentClass = candidateKey.front();
        auto result = pathwayStorage.states.try_emplace(
            std::move(candidateKey),
            pathwayTransposition{sumDupBonds, nullptr}
        );
        inserted = result.second;
        bestSumDupBonds = &result.first->second.bestSumDupBonds;

        if (inserted)
        {
            try
            {
                existingPath = pathwayStorage.createPath(
                    retainedFragmentClass,
                    parent.apPtr,
                    bitsetHashTable.find(matching.first)->second.second,
                    bitsetHashTable.find(matching.second)->second.second
                );
            }
            catch (...)
            {
                pathwayStorage.states.erase(result.first);
                throw;
            }
            result.first->second.path = existingPath;
            candidate.apPtr = existingPath;
        }
        else existingPath = result.first->second.path;
    }

#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (searchTelemetryEnabled) [[unlikely]]
    {
        if (inserted) ++searchTelemetry.counters.assemblyCacheMisses;
        else ++searchTelemetry.counters.assemblyCacheHits;
    }
#endif

    if (!inserted)
    {
        if (sumDupBonds <= *bestSumDupBonds)
        {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
            if (searchTelemetryEnabled) [[unlikely]]
                ++searchTelemetry.counters.assemblyCachePrunedHits;
#endif
            return true;
        }

#ifdef ASSEMBLY_ENABLE_TELEMETRY
        if (searchTelemetryEnabled) [[unlikely]]
            ++searchTelemetry.counters.assemblyCacheUpdatedHits;
#endif
        *bestSumDupBonds = sumDupBonds;
        if (searchStorage.tracksPath())
        {
            candidate.apPtr = existingPath;
            setAssemblyPathStep(*existingPath, parent.apPtr, matching);
        }
    }

    dagRecursiveAssemblyWithWorkspaceImpl<equivalenceMode>(
        candidate,
        AI,
        fragmentationWorkspace,
        searchStorage
    );
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
template<matchingEquivalenceMode equivalenceMode>
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
    recordImprovedAssemblyIndex(input, AI, searchStorage);
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
    candidate.reserveFragments(input.masks.size() + 2);
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        if (searchShouldStop()) return;
        dagDuplicateClassLevel &stmap = stmapVector[j];
        vector<EdgeMask> stmapMaskList(input.masks.size(), 0);
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
                                            input.masks.size() * input.masks.size(),
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
                                        firstFragment * input.masks.size() +
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
                        vi candidateKey;
                        if (!canoniseAssemblyStateAndBuildKey(
                            candidate,
                            candidateKey
                        )) return false;
                        if (!continueAssemblySearchWithWorkspace<
                            equivalenceMode
                        >(
                            input,
                            candidate,
                            matching,
                            std::move(candidateKey),
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

void dagRecursiveAssemblyWithWorkspace(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    if (!fragmentationWorkspace.homogeneousPathEdgePositions.empty())
    {
        dagRecursiveAssemblyWithWorkspaceImpl<
            matchingEquivalenceMode::homogeneousPath
        >(input, AI, fragmentationWorkspace, searchStorage);
    }
    else
    {
        dagRecursiveAssemblyWithWorkspaceImpl<
            matchingEquivalenceMode::none
        >(input, AI, fragmentationWorkspace, searchStorage);
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
template<matchingEquivalenceMode equivalenceMode>
void initialRecursiveAssemblyWithWorkspaceImpl(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    assemblySearchStorage &searchStorage
)
{
    recordImprovedAssemblyIndex(input, AI, searchStorage);
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
    candidate.reserveFragments(input.masks.size() + 2);
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
                candidate.clearFragments();
                fragmentAssemblyStateWithoutCanonisationWithWorkspace(
                    input,
                    matching,
                    candidate,
                    fragmentationWorkspace
                );
                if (searchShouldStop()) return false;
                int sumDupBonds =
                    input.sumDupBonds + matching.maxFragSize - 1;
                candidate.sumDupBonds = sumDupBonds;
                if (candidate.lowBoundAI() < AI)
                {
                    vi candidateKey;
                    if (!canoniseAssemblyStateAndBuildKey(
                        candidate,
                        candidateKey
                    )) return false;
                    if (!continueAssemblySearchWithWorkspace<
                        equivalenceMode
                    >(
                        input,
                        candidate,
                        matching,
                        std::move(candidateKey),
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
            matchingEquivalenceMode::homogeneousPath
        >(input, AI, fragmentationWorkspace, searchStorage);
    }
    else
    {
        initialRecursiveAssemblyWithWorkspaceImpl<
            matchingEquivalenceMode::none
        >(input, AI, fragmentationWorkspace, searchStorage);
    }
}

bool runImprovedAssemblySearch(
    assemblyState &root,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace,
    vector<edgeL> &removedEdges,
    ofstream &ofs,
    assemblySearchStorage &searchStorage
)
{
    vi rootKey = root.assemblyHashCalculator();
    if (searchStorage.tracksPath())
    {
        pathwayAssemblySearchStorage &pathwayStorage = *searchStorage.pathway;
        assemblyPath *rootPath = pathwayStorage.createPath(
            rootKey.front(),
            nullptr,
            0,
            0
        );
        pathwayStorage.states.emplace(
            std::move(rootKey),
            pathwayTransposition{0, rootPath}
        );
        root.apPtr = rootPath;
    }
    else
    {
        searchStorage.scoreOnly->states.emplace(std::move(rootKey), 0);
    }

    initialRecursiveAssemblyWithWorkspace(
        root,
        AI,
        fragmentationWorkspace,
        searchStorage
    );

#ifdef ASSEMBLY_ENABLE_TELEMETRY
    setSearchTelemetryPhase(SearchTelemetryPhase::output);
#endif
    // Keep the primary result line machine-readable even for partial searches.
    ofs << compensateDisjointAssemblyIndex(AI) << '\n';
    if (runtimeLimitReached)
    {
        cout << "status: runtime limit reached\n";
        ofs << "status: runtime limit reached\n";
    }
    if (enumerationLimitReached)
    {
        cout << "status: enumeration limit reached\n";
        ofs << "status: enumeration limit reached\n";
    }

    if (searchStorage.tracksPath())
    {
        const bool recovered = recoverPathway2(removedEdges);
        minAssemblyPath = nullptr;
        return recovered;
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
    minAssemblyPath = nullptr;
    bitsetHashTable.clear();
    graphHashMap.clear();
    clearTreeCanonInterner();
    DAG.clear();
    intermediateMAs.clear();
    totalBonds = mg.totalBonds;
    originalEdgeList = mg.writeEdgeList();
    disjointFragments = mg.disjointFragments();
    originalMolecule = mg;
    vector<edgeL> removedEdges;
    targetMolecule = preprocessWriteback(mg, removedEdges);
    univEdgeList = targetMolecule.writeEdgeList();
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
    as.appendFragment(allEdges, static_cast<int>(univEdgeList.size()));
    int AI = std::numeric_limits<int>::max();
    if (isPathway)
    {
        pathwayAssemblySearchStorage pathwayStorage;
        assemblySearchStorage searchStorage{nullptr, &pathwayStorage};
        return runImprovedAssemblySearch(
            as,
            AI,
            fragmentationWorkspace,
            removedEdges,
            ofs,
            searchStorage
        );
    }
    scoreOnlyAssemblySearchStorage scoreOnlyStorage;
    assemblySearchStorage searchStorage{&scoreOnlyStorage, nullptr};
    return runImprovedAssemblySearch(
        as,
        AI,
        fragmentationWorkspace,
        removedEdges,
        ofs,
        searchStorage
    );
}
