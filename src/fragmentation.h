#pragma once

/**
 * @brief Splits an assembly state without canonicalising the resulting masks
 *
 * @param target The assembly state to be fragmented
 * @param matching The duplicate pair. matching.first is retained and matching.second is deleted
 * @param duplicateCanonicalId Canonical ID of the enclosing duplicate class
 * @param result The resulting assembly state
 * @param workspace Buffers reused by successive fragmentation calls
 *
 * The caller may evaluate bitset-only bounds on the raw result. Unknown IDs
 * must be resolved before the state is hashed or recursively enumerated.
 */
void fragmentAssemblyStateWithoutCanonisationWithWorkspace(
    assemblyState &target,
    validMatchings &matching,
    int duplicateCanonicalId,
    assemblyState &result,
    ufdsMaskWorkspace &workspace
)
{
    workspace.beginFragmentation();
    vector<assemblyFragment> &fragments = target.fragments;
    EdgeMask f1 = matching.first, f2 = matching.second;
    const bool same =
        matching.firstFragmentIndex == matching.secondFragmentIndex;
    result.appendFragment(
        f1,
        matching.maximumFragmentSize,
        duplicateCanonicalId,
        true
    );
    if (same)
    {
        EdgeMask resultMask = fragments[matching.firstFragmentIndex].mask;
        resultMask ^= f1;
        resultMask ^= f2;
        ufdsMaskConstructWithWorkspace(
            resultMask,
            result.fragments,
            workspace
        );
    }
    else
    {
        EdgeMask resultMask1 = fragments[matching.firstFragmentIndex].mask;
        resultMask1 ^= f1;
        ufdsMaskConstructWithWorkspace(
            resultMask1,
            result.fragments,
            workspace
        );
        EdgeMask resultMask2 = fragments[matching.secondFragmentIndex].mask;
        resultMask2 ^= f2;
        ufdsMaskConstructWithWorkspace(
            resultMask2,
            result.fragments,
            workspace
        );
    }
    for (size_t fragmentIndex = 0;
         fragmentIndex < fragments.size();
         fragmentIndex++)
    {
        const assemblyFragment &fragment = fragments[fragmentIndex];
        if (
            fragmentIndex !=
                static_cast<size_t>(matching.firstFragmentIndex) &&
            fragmentIndex !=
                static_cast<size_t>(matching.secondFragmentIndex) &&
            fragment.mask != 0
        )
        {
            // DAG enumeration may trim a parent mask enough to disconnect it.
            // Preserve proved metadata; uncertain masks must be split.
            if (!fragment.connected)
            {
                ufdsMaskConstructWithWorkspace(
                    fragment.mask,
                    result.fragments,
                    workspace
                );
            }
            else result.appendFragment(fragment);
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
    IntegerVector &key,
    ufdsMaskWorkspace &workspace
)
{
    if (searchShouldStop()) return false;
    key.resize(target.fragments.size());
    for (
        size_t fragmentIndex = 0;
        fragmentIndex < target.fragments.size();
        fragmentIndex++
    )
    {
        assemblyFragment &fragment = target.fragments[fragmentIndex];
        if (fragment.canonicalId < 0)
        {
            fragment.canonicalId = canonise(fragment.mask);
            if (searchShouldStop()) return false;
        }
        key[fragmentIndex] = fragment.canonicalId;
    }
    if (key.size() > 1) sort(key.begin() + 1, key.end());
    workspace.cacheCanonicalIds(target.fragments);
    return true;
}
