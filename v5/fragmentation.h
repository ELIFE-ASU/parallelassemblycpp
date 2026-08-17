/**
 * @brief Splits an assembly state without canonicalising the resulting masks
 *
 * @param _target The assembly state to be fragmented
 * @param matching The duplicate pair. matching.first is retained and matching.second is deleted
 * @param _result The resulting assembly state
 * @param workspace Buffers reused by successive fragmentation calls
 *
 * The caller may evaluate bitset-only bounds on the raw result. Every mask must
 * be canonicalised before the state is hashed or recursively enumerated.
 */
void fragmentAssemblyStateWithoutCanonisationWithWorkspace(
    assemblyState &_target,
    validMatchings &matching,
    assemblyState &_result,
    ufdsMaskWorkspace &workspace
)
{
    vector<EdgeMask> &masks = _target.masks;
    EdgeMask f1 = matching.first, f2 = matching.second;
    const bool same = matching.frag1 == matching.frag2;
    _result.masks.push_back(f1);
    if (same)
    {
        EdgeMask resultMask = masks[matching.frag1];
        resultMask ^= f1;
        resultMask ^= f2;
        ufdsMaskConstructWithWorkspace(resultMask, _result.masks, workspace);
    }
    else
    {
        EdgeMask resultMask1 = masks[matching.frag1];
        resultMask1 ^= f1;
        ufdsMaskConstructWithWorkspace(resultMask1, _result.masks, workspace);
        EdgeMask resultMask2 = masks[matching.frag2];
        resultMask2 ^= f2;
        ufdsMaskConstructWithWorkspace(resultMask2, _result.masks, workspace);
    }
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (
            i != static_cast<size_t>(matching.frag1) &&
            i != static_cast<size_t>(matching.frag2) &&
            masks[i] != 0
        )
        {
            // DAG enumeration may trim a parent mask enough to disconnect it.
            // Cached masks are known connected; uncached masks must be split.
            if (bitsetHashTable.count(masks[i]) == 0)
            {
                ufdsMaskConstructWithWorkspace(
                    masks[i],
                    _result.masks,
                    workspace
                );
            }
            else _result.masks.push_back(masks[i]);
        }
    }
}

/**
 * @brief Canonicalise a bound-surviving state and build its transposition key
 *
 * @param target The raw fragmented state
 * @param key Canonical fragment indices, with the retained fragment first
 * @return false when the search should stop, otherwise true
 */
bool canoniseAssemblyStateAndBuildKey(assemblyState &target, vi &key)
{
    if (searchShouldStop()) return false;
    key.resize(target.masks.size());
    for (size_t i = 0; i < target.masks.size(); i++)
    {
        key[i] = canonise(target.masks[i]);
        if (searchShouldStop()) return false;
    }
    if (key.size() > 1) sort(key.begin() + 1, key.end());
    return true;
}

/**
 * @brief Empty the hash table
 * 
 */
void clearPathMap()
{
    for (const apWrapper &wrapper : pathAssemblyMap)
    {
        delete wrapper.ap;
    }
    pathAssemblyMap.clear();
}
