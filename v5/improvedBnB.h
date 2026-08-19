/**
 * @brief Enumerate all subgraphs during the initial phase of the pathway algorithm. See Seet et al section 4.3 Duplicate Enumeration
 *
 * @param _target The initial assembly state
 * @param stmapVector The matchings found
 * @return true if any matchings found
 * @return false if no matchings found
 */
bool initialRecursiveEnumeration(assemblyState &_target, vector<map<int, initialDuplicateSet> > &stmapVector)
{
    vector<initialDagLevel> tempDag(2);
    vector<EdgeMask> &masks = _target.masks;
    bool alive = 0;
    size_t currSize = 1;
    vector<initialPotentialDuplicate> prevML;
    size_t retainedStateCount = 0;

    // Retain the one-edge DAG states first so they count toward ENUM_MAX too.
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStop()) return false;
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
            }
        }
    }

    // Generate the first multi-edge frontier only after all base states exist.
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStop()) return false;
            if (masks[i][j] == 0) continue;

            initialPotentialDuplicate m(j, masks[i], i);
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
        map<int, initialDuplicateSet> &stmap = stmapVector.back();
        active = 0;
        vector<initialPotentialDuplicate> currML;
        for (size_t i = 0; i < prevML.size(); i++)
        {
            if (searchShouldStop()) return false;
            initialPotentialDuplicate &m = prevML[i];

            int s = canonise(m.mask);
            if (searchShouldStop()) return false;
            auto entry = stmap.try_emplace(
                s,
                currSize + 1,
                masks.size()
            ).first;
            entry->second.insert(m);
        }
        tempDag.resize(tempDag.size() + 1);
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (searchShouldStop()) return false;
            initialDuplicateSet &ss = it->second;
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
int dagRecursiveEnumeration(assemblyState &_target, vector<map<int, dagDuplicateSet> > &stmapVector,
    vector<vector<EdgeMask> > &targetMasks)
{
    if (searchShouldStop()) return 0;
    int ordinal = std::numeric_limits<int>::max();
    const auto ordinalEntry = bitsetHashTable.find(_target.masks.front());
    if (ordinalEntry != bitsetHashTable.end()) ordinal = ordinalEntry->second.first;
    vector<EdgeMask> &masks = _target.masks;
    size_t currSize = 1;
    
    stmapVector.resize(1);
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStop()) return 0;
            if (masks[i][j] != 0)
            {
                EdgeMask b = 0; b.set(j);
                potentialDuplicate m(b, i, j);
                dagGenerate(m, stmapVector[0], masks[i], currSize, ordinal, masks.size());
                if (searchShouldStop()) return 0;
            }
        }
    }
    bool active = 1, overweight = 0, last = 0;
    while (active)
    {
        if (searchShouldStop()) return 0;
        vector<EdgeMask> targetMask(masks.size(), 0);
        active = 0;
        stmapVector.emplace_back();
        map<int, dagDuplicateSet> &stmap = stmapVector[stmapVector.size() - 2];
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (searchShouldStop()) return 0;
            dagDuplicateSet &ss = it->second;
            if (ss.isValid())
            {
                if (searchShouldStop()) return 0;
                active |= dagDuplicateGenerator(ss, stmapVector.back(), targetMask, masks, ordinal, overweight, last);
                if (searchShouldStop()) return 0;
            }
        }
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

void dagRecursiveAssemblyWithWorkspace(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace
);

void recordImprovedAssemblyIndex(assemblyState &input, int &AI)
{
    const int candidate = input.AI();
    if (candidate >= AI) return;

    AI = candidate;
    minAssemblyPath = input.apPtr;
    const unsigned long long time = elapsedClockTicks();
    cout << "time: " << time << " min AI found so far: " << AI << '\n';
    if (writeIntermediateMAs) intermediateMAs.emplace_back(time, AI);
}

static std::unique_ptr<assemblyPath> createAssemblyPath(
    vi key,
    assemblyPath *parent,
    validMatchings &matching,
    int sumDupBonds
)
{
    auto path = std::make_unique<assemblyPath>();
    path->key = std::move(key);
    path->parent = parent;
    path->sumDupBonds = sumDupBonds;
    path->match = bitsetHashTable[matching.first].second;
    path->duplicate = bitsetHashTable[matching.second].second;
    return path;
}

bool continueAssemblySearchWithWorkspace(
    assemblyState &parent,
    assemblyState &candidate,
    validMatchings &matching,
    vi candidateKey,
    int sumDupBonds,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace
)
{
    // For larger keys, one lookup plus a miss-only copy beats hashing twice.
    constexpr size_t singleInsertMinFragments = 5;
    assemblyPath *existingPath = nullptr;
    {
        assemblyPath probe{};
        probe.key = std::move(candidateKey);
        if (searchShouldStop()) return false;

        if (probe.key.size() < singleInsertMinFragments)
        {
            auto entry = pathAssemblyMap.find(apWrapper{&probe});
            if (entry == pathAssemblyMap.end())
            {
                auto path = createAssemblyPath(
                    std::move(probe.key),
                    parent.apPtr,
                    matching,
                    sumDupBonds
                );
                auto [position, inserted] =
                    pathAssemblyMap.insert(apWrapper{path.get()});
                if (inserted)
                {
                    candidate.apPtr = path.get();
                    path.release();
                }
                else existingPath = position->ap;
            }
            else existingPath = entry->ap;
        }
        else
        {
            auto [entry, inserted] = pathAssemblyMap.insert(apWrapper{&probe});
            if (inserted)
            {
                std::unique_ptr<assemblyPath> path;
                try
                {
                    path = createAssemblyPath(
                        probe.key,
                        parent.apPtr,
                        matching,
                        sumDupBonds
                    );
                }
                catch (...)
                {
                    pathAssemblyMap.erase(entry);
                    throw;
                }

                // The copied key keeps the set's hash and equality unchanged.
                entry->ap = path.get();
                candidate.apPtr = path.get();
                path.release();
            }
            else existingPath = entry->ap;
        }

        if (existingPath != nullptr)
        {
            if (sumDupBonds <= existingPath->sumDupBonds) return true;

            candidate.apPtr = existingPath;
            existingPath->sumDupBonds = sumDupBonds;
            existingPath->match = bitsetHashTable[matching.first].second;
            existingPath->duplicate = bitsetHashTable[matching.second].second;
            existingPath->parent = parent.apPtr;
        }
    }

    dagRecursiveAssemblyWithWorkspace(
        candidate,
        AI,
        fragmentationWorkspace
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
void dagRecursiveAssemblyWithWorkspace(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace
)
{
    // Count-only pair pruning pays for itself once residual splitting is no
    // longer dominated by the fixed cost of the bound lookup.
    constexpr size_t pairBoundMinimumMoleculeEdges = 27;
    const bool usePairBound =
        fragmentationWorkspace.edgeCount >= pairBoundMinimumMoleculeEdges;
    recordImprovedAssemblyIndex(input, AI);
    if (searchShouldStop()) return;

    vector<map<int, dagDuplicateSet> > stmapVector;
    vector<vector<EdgeMask> > targetMasks;
    int maxFragSize = dagRecursiveEnumeration(input, stmapVector, targetMasks);

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
        map<int, dagDuplicateSet> &stmap = stmapVector[j];
        vector<EdgeMask> stmapMaskList(input.masks.size(), 0);
        vi unrestrictedParentTotals;
        vi pairGenericBoundCache;
        EdgeMask maskM = 0;
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (searchShouldStop()) return;
            dagDuplicateSet &ss = it->second;
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
                const bool completed = ss.visitMatchingsInReverse(
                [&](validMatchings &matching)
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

                    candidate.clearFragments();
                    fragmentAssemblyStateWithoutCanonisationWithWorkspace(
                        input,
                        matching,
                        candidate,
                        fragmentationWorkspace
                    );
                    if (searchShouldStop()) return false;

                    int sumDupBonds = input.sumDupBonds + matching.maxFragSize - 1;
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
                        if (!canoniseAssemblyStateAndBuildKey(candidate, candidateKey))
                            return false;
                        if (!continueAssemblySearchWithWorkspace(
                            input,
                            candidate,
                            matching,
                            std::move(candidateKey),
                            sumDupBonds,
                            AI,
                            fragmentationWorkspace
                        )) return false;
                    }
                    return true;
                });
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
void initialRecursiveAssemblyWithWorkspace(
    assemblyState &input,
    int &AI,
    ufdsMaskWorkspace &fragmentationWorkspace
)
{
    recordImprovedAssemblyIndex(input, AI);
    if (searchShouldStop()) return;

    vector<map<int, initialDuplicateSet> > stmapVector;
    const bool hasInitialMatchings = initialRecursiveEnumeration(input, stmapVector);
    if (enumerationLimitReached || searchShouldStop()) return;

    if (!hasInitialMatchings) return;
    assemblyState candidate;
    candidate.reserveFragments(input.masks.size() + 2);
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        if (searchShouldStop()) return;
        map<int, initialDuplicateSet> &stmap = stmapVector[j];
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (searchShouldStop()) return;
            initialDuplicateSet &ss = it->second;
            const bool completed = ss.visitMatchingsInReverse(
            [&](validMatchings &matching)
            {
                candidate.clearFragments();
                fragmentAssemblyStateWithoutCanonisationWithWorkspace(
                    input,
                    matching,
                    candidate,
                    fragmentationWorkspace
                );
                if (searchShouldStop()) return false;
                int sumDupBonds = input.sumDupBonds + matching.maxFragSize - 1;
                candidate.sumDupBonds = sumDupBonds;
                if (candidate.lowBoundAI() < AI)
                {
                    vi candidateKey;
                    if (!canoniseAssemblyStateAndBuildKey(candidate, candidateKey))
                        return false;
                    if (!continueAssemblySearchWithWorkspace(
                        input,
                        candidate,
                        matching,
                        std::move(candidateKey),
                        sumDupBonds,
                        AI,
                        fragmentationWorkspace
                    )) return false;
                }
                return true;
            });
            if (!completed) return;
        }
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
    runtimeLimitReached = false;
    enumerationLimitReached = false;
    clearPathMap();
    minAssemblyPath = nullptr;
    bitsetHashTable.clear();
    graphHashMap.clear();
    DAG.clear();
    intermediateMAs.clear();
    totalBonds = mg.totalBonds;
    originalEdgeList = mg.writeEdgeList();
    disjointFragments = mg.disjointFragments();
    originalMolecule = mg;
    vector<edgeL> removedEdges;
    targetMolecule = preprocessWriteback(mg, removedEdges);
    univEdgeList = targetMolecule.writeEdgeList();
    EdgeMask::configure(univEdgeList.size());
    AtomMask::configure(targetMolecule.mg.size());
    ufdsMaskWorkspace fragmentationWorkspace(
        targetMolecule.mg.size(),
        univEdgeList.size()
    );
    allEdges.reset();
    for (size_t i = 0; i < univEdgeList.size(); i++) allEdges.set(i);
    assemblyState as;
    as.appendFragment(allEdges, static_cast<int>(univEdgeList.size()));
    auto rootPath = std::make_unique<assemblyPath>();
    apWrapper ap{rootPath.get()};
    ap.ap->parent = nullptr;
    ap.ap->key = as.assemblyHashCalculator();
    ap.ap->sumDupBonds = 0;
    pathAssemblyMap.insert(ap);
    as.apPtr = rootPath.release();
    int AI = std::numeric_limits<int>::max();
    initialRecursiveAssemblyWithWorkspace(as, AI, fragmentationWorkspace);

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
    if (isPathway) return recoverPathway2(removedEdges);
    return true;
}
