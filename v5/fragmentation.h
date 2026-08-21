/**
 * @brief Splits an assembly state without canonicalising the resulting masks
 *
 * @param _target The assembly state to be fragmented
 * @param matching The duplicate pair. matching.first is retained and matching.second is deleted
 * @param duplicateCanonicalId Canonical ID of the enclosing duplicate class
 * @param _result The resulting assembly state
 * @param workspace Buffers reused by successive fragmentation calls
 *
 * The caller may evaluate bitset-only bounds on the raw result. Unknown IDs
 * must be resolved before the state is hashed or recursively enumerated.
 */
void fragmentAssemblyStateWithoutCanonisationWithWorkspace(
    assemblyState &_target,
    validMatchings &matching,
    int duplicateCanonicalId,
    assemblyState &_result,
    ufdsMaskWorkspace &workspace
)
{
    workspace.beginFragmentation();
    vector<assemblyFragment> &fragments = _target.fragments;
    EdgeMask f1 = matching.first, f2 = matching.second;
    const bool same = matching.frag1 == matching.frag2;
    _result.appendFragment(
        f1,
        matching.maxFragSize,
        duplicateCanonicalId,
        true
    );
    if (same)
    {
        EdgeMask resultMask = fragments[matching.frag1].mask;
        resultMask ^= f1;
        resultMask ^= f2;
        ufdsMaskConstructWithWorkspace(
            resultMask,
            _result.fragments,
            workspace
        );
    }
    else
    {
        EdgeMask resultMask1 = fragments[matching.frag1].mask;
        resultMask1 ^= f1;
        ufdsMaskConstructWithWorkspace(
            resultMask1,
            _result.fragments,
            workspace
        );
        EdgeMask resultMask2 = fragments[matching.frag2].mask;
        resultMask2 ^= f2;
        ufdsMaskConstructWithWorkspace(
            resultMask2,
            _result.fragments,
            workspace
        );
    }
    for (size_t i = 0; i < fragments.size(); i++)
    {
        const assemblyFragment &fragment = fragments[i];
        if (
            i != static_cast<size_t>(matching.frag1) &&
            i != static_cast<size_t>(matching.frag2) &&
            fragment.mask != 0
        )
        {
            // DAG enumeration may trim a parent mask enough to disconnect it.
            // Preserve proved metadata; uncertain masks must be split.
            if (!fragment.connected)
            {
                ufdsMaskConstructWithWorkspace(
                    fragment.mask,
                    _result.fragments,
                    workspace
                );
            }
            else _result.appendFragment(fragment);
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
bool canoniseAssemblyStateAndBuildKey(
    assemblyState &target,
    vi &key,
    ufdsMaskWorkspace &workspace
)
{
    if (searchShouldStop()) return false;
    key.resize(target.fragments.size());
    for (size_t i = 0; i < target.fragments.size(); i++)
    {
        assemblyFragment &fragment = target.fragments[i];
        if (fragment.canonicalId < 0)
        {
            fragment.canonicalId = canonise(fragment.mask);
            if (searchShouldStop()) return false;
        }
        key[i] = fragment.canonicalId;
    }
    if (key.size() > 1) sort(key.begin() + 1, key.end());
    workspace.cacheCanonicalIds(target.fragments);
    return true;
}
