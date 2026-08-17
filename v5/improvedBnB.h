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
    vector<std::unordered_map<standardBitset, pair<int, vector<standardBitset> > > > tempDag(2);
    vector<standardBitset> &masks = _target.masks;
    bool alive = 0;
    size_t currSize = 1;
    vector<initialPotentialDuplicate> prevML;
    std::unordered_set<standardBitset> maskMap;

    // Retain the one-edge DAG states first so they count toward ENUM_MAX too.
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStop()) return false;
            if (masks[i][j] != 0)
            {
                pair<int, vector<standardBitset> > p; p.first = -1;
                standardBitset b = 0; b.set(j);
                if (!reserveInitialDagMask(maskMap, b))
                {
                    return false;
                }
                tempDag[0][b] = p;
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
            if (!m.generateDAG(prevML, maskMap, tempDag))
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
                alive |= ss.dagPopulator(currML, maskMap, tempDag);
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
    vector<vector<standardBitset> > &targetMasks)
{
    if (searchShouldStop()) return 0;
    int ordinal = std::numeric_limits<int>::max();
    const auto ordinalEntry = bitsetHashTable.find(_target.masks.front());
    if (ordinalEntry != bitsetHashTable.end()) ordinal = ordinalEntry->second.first;
    vector<standardBitset> &masks = _target.masks;
    size_t currSize = 1;
    
    stmapVector.resize(1);
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (searchShouldStop()) return 0;
            if (masks[i][j] != 0)
            {
                standardBitset b = 0; b.set(j);
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
        vector<standardBitset> targetMask(masks.size(), 0);
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
    }
    return currSize;
}


/**
 * @brief This function returns the minimum possible sum of duplicate bonds that can be found if only fragments equal
 * to the size of the largest duplicatable subgraph are considered. Used to tighten the bound after fragmentation
 * @param target The target assembly state
 * @param matchMask The bitset of all graphs isomorphic to the matching
 * @param maxFragMask The bitset of all duplicatable subgraphs with the same bitset count as the matching
 * @return int the sum of duplicate bonds calculated by this cutoff function
 */
int postFragmentationCutoff(assemblyState &target, standardBitset & matchMask, standardBitset &maxFragMask)
{
    if (target.masks.size() < 2) return 0;
    int maxFragSize = target.masks[0].count();

    /// mainSizeList is the vector of integers corresponding to the sizes of each fragment
    /// matchSizeList is the vector of integers corresponding to the sizes of each matching fragment bitset
    /// maxFragSizeList is the vector of integers corrsponding to the sizes of each maximal size fragment bitset
    vi mainSizeList, matchSizeList, maxFragSizeList;

    mainSizeList.resize(target.masks.size());
    matchSizeList.resize(target.masks.size());
    maxFragSizeList.resize(target.masks.size());

    mainSizeList[0] = maxFragSize;
    matchSizeList[0] = maxFragSize;
    maxFragSizeList[0] = maxFragSize;
    for (size_t i = 1; i < target.masks.size(); i++)
    {
        mainSizeList[i] = target.masks[i].count();
        matchSizeList[i] = (target.masks[i] & matchMask).count();
        maxFragSizeList[i] = (target.masks[i] & maxFragMask).count();
    }
    int matchDB = target.maxDupBonds(mainSizeList, maxFragSize, matchSizeList);
    int maxFragDB = target.maxDupBonds(mainSizeList, maxFragSize, maxFragSizeList) - 1;
    return max(matchDB, maxFragDB);
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
        probe.key = candidate.assemblyHashCalculator();
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
    recordImprovedAssemblyIndex(input, AI);
    if (searchShouldStop()) return;

    vector<map<int, dagDuplicateSet> > stmapVector;
    vector<vector<standardBitset> > targetMasks;
    int maxFragSize = dagRecursiveEnumeration(input, stmapVector, targetMasks);

    if (searchShouldStop()) return;

    /// Find the fragment-size-specific AI lower bounds
    vi fragSizeList, fragSizeListMax, sizeList(input.masks.size());
    input.maxDupBonds(fragSizeList, maxFragSize, targetMasks);
    if (searchShouldStop()) return;

    /// Establish the max duplicate cutoff based on the maximum fragment size
    fragSizeListMax.resize(fragSizeList.size());
    fragSizeListMax[0] = fragSizeList[0];
    for (size_t i = 1; i < fragSizeList.size(); i++)
    {
        if (searchShouldStop()) return;
        if (fragSizeList[i] > fragSizeListMax[i - 1])
            fragSizeListMax[i] = fragSizeList[i];
        else fragSizeListMax[i] = fragSizeListMax[i - 1];
    }
    for (size_t i = 0; i < input.masks.size(); i++)
    {
        if (searchShouldStop()) return;
        sizeList[i] = input.masks[i].count();
    }

    /// Begin iterating through the enumerated duplicatable fragments
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        if (searchShouldStop()) return;
        map<int, dagDuplicateSet> &stmap = stmapVector[j];
        vector<standardBitset> stmapMaskList(input.masks.size(), 0);
        standardBitset maskM = 0;
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (searchShouldStop()) return;
            dagDuplicateSet &ss = it->second;
            if (!ss.dead)
            {
            vector<validMatchings> matchings;
            standardBitset maskC = 0;
            
            for (size_t i = 0; i < ss.maskList.size(); i++)
            {
                if (searchShouldStop()) return;
                maskC |= ss.maskList[i];
                stmapMaskList[i] |= ss.maskList[i];
            }
            maskM |= maskC;
            int dupBondsMaxFrag = input.maxDupBonds(sizeList, ss.size, stmapMaskList);

            int temp = fragSizeListMax[ss.size - 2] - 1;

            if (ss.size > 2) temp = max(temp, fragSizeListMax[ss.size - 3]);

            /// Initial branch-and-bound before the fragmentation step
            int earlySDP = max(dupBondsMaxFrag, temp);

            int earlyAIBound = totalBonds - input.sumDupBonds - 1 - earlySDP;
            if (earlyAIBound < AI)
            {
                if (ss.list.size() > 1)
                {
                    ss.generateMatchings(matchings);
                    if (searchShouldStop()) return;
                }
                for (int i = matchings.size() - 1; i >= 0; i--)
                {
                    if (searchShouldStop()) return;
                    assemblyState as;
                    fragmentAssemblyStateWithWorkspace(
                        input,
                        matchings[i],
                        as,
                        fragmentationWorkspace
                    );
                    if (searchShouldStop()) return;

                    int sumDupBonds = input.sumDupBonds + matchings[i].maxFragSize - 1;
                    as.sumDupBonds = sumDupBonds;
                    int fragmentationCutoff = postFragmentationCutoff(as, maskC, maskM);
                    if (searchShouldStop()) return;
                    if (as.lowBoundAI(matchings[i].maxFragSize, fragmentationCutoff) < AI)
                    {
                        if (!continueAssemblySearchWithWorkspace(
                            input,
                            as,
                            matchings[i],
                            sumDupBonds,
                            AI,
                            fragmentationWorkspace
                        )) return;
                    }
                }
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
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        if (searchShouldStop()) return;
        map<int, initialDuplicateSet> &stmap = stmapVector[j];
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (searchShouldStop()) return;
            initialDuplicateSet &ss = it->second;
            vector<validMatchings> matchings;
            if (ss.list.size() > 1)
            {
                ss.generateMatchings(matchings);
                if (searchShouldStop()) return;
            }
            for (int i = matchings.size() - 1; i >= 0; i--)
            {
                if (searchShouldStop()) return;
                assemblyState as;
                fragmentAssemblyStateWithWorkspace(
                    input,
                    matchings[i],
                    as,
                    fragmentationWorkspace
                );
                if (searchShouldStop()) return;
                int sumDupBonds = input.sumDupBonds + matchings[i].maxFragSize - 1;
                as.sumDupBonds = sumDupBonds;
                if (as.lowBoundAI() < AI)
                {
                    if (!continueAssemblySearchWithWorkspace(
                        input,
                        as,
                        matchings[i],
                        sumDupBonds,
                        AI,
                        fragmentationWorkspace
                    )) return;
                }
            }
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
    ufdsMaskWorkspace fragmentationWorkspace(
        targetMolecule.mg.size(),
        univEdgeList.size()
    );
    allEdges = 0;
    for (size_t i = 0; i < univEdgeList.size(); i++) allEdges.set(i);
    assemblyState as; as.masks.push_back(allEdges);
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
