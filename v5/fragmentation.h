/**
 * @brief Splits the assembly state into fragments after a given duplicate is removed using the disjoint-set data structure
 *
 * @param _target The assembly state to be fragmented
 * @param matching The duplicate pair. matching.first is retained and matching.second is deleted
 * @param _result The resulting assembly state
 * @param workspace Buffers reused by successive fragmentation calls
 */
void fragmentAssemblyStateWithWorkspace(
    assemblyState &_target,
    validMatchings &matching,
    assemblyState &_result,
    ufdsMaskWorkspace &workspace
)
{
    vector<standardBitset> &masks = _target.masks;
    standardBitset f1 = matching.first, f2 = matching.second;
    const bool same = matching.frag1 == matching.frag2;
    _result.masks.push_back(f1);
    if (same)
    {
        standardBitset resultMask = masks[matching.frag1];
        resultMask ^= f1;
        resultMask ^= f2;
        ufdsMaskConstructWithWorkspace(resultMask, _result.masks, workspace);
    }
    else
    {
        standardBitset resultMask1 = masks[matching.frag1];
        resultMask1 ^= f1;
        ufdsMaskConstructWithWorkspace(resultMask1, _result.masks, workspace);
        standardBitset resultMask2 = masks[matching.frag2];
        resultMask2 ^= f2;
        ufdsMaskConstructWithWorkspace(resultMask2, _result.masks, workspace);
    }
    for (size_t i = 0; i < _result.masks.size(); i++)
    {
        canonise(_result.masks[i]);
    }
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (
            i != static_cast<size_t>(matching.frag1) &&
            i != static_cast<size_t>(matching.frag2) &&
            masks[i] != 0
        )
        {
            if (bitsetHashTable.count(masks[i]) == 0)
            {
                const size_t firstNewMask = _result.masks.size();
                ufdsMaskConstructWithWorkspace(
                    masks[i],
                    _result.masks,
                    workspace
                );
                for (size_t j = firstNewMask; j < _result.masks.size(); j++)
                {
                    canonise(_result.masks[j]);
                }
            }
            else _result.masks.push_back(masks[i]);
        }
    }
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
