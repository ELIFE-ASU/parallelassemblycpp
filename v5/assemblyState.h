/**
 * @brief Assembly state data structure. Records the current state of this assembly pathway
 */
struct assemblyState
{   
    /// @brief each mask represents a separate fragment as a boolean edge list
    vector<EdgeMask> masks;
    /// @brief edge count parallel to each entry in masks
    vi edgeCounts;
    /// @brief number of duplicated bonds
    int sumDupBonds = 0;
    void appendFragment(const EdgeMask &mask, int edgeCount)
    {
        masks.push_back(mask);
        edgeCounts.push_back(edgeCount);
    }

    void reserveFragments(size_t capacity)
    {
        masks.reserve(capacity);
        edgeCounts.reserve(capacity);
    }

    void clearFragments()
    {
        masks.clear();
        edgeCounts.clear();
        sumDupBonds = 0;
    }

    static int fixedSizeDupBondsForFragment(
        int fragmentEdges,
        int duplicateSize,
        int eligibleEdges
    )
    {
        const int completeGroups = eligibleEdges / duplicateSize;
        const int adjustedSize = completeGroups * duplicateSize;
        const int remainingSize = fragmentEdges - adjustedSize;
        int result = completeGroups * (duplicateSize - 1);
        result += remainingSize - remainingSize / (duplicateSize - 1);
        if (remainingSize % (duplicateSize - 1) != 0) result--;
        return result;
    }

    static int unrestrictedDupBondsForFragment(
        int fragmentEdges,
        int duplicateSize
    )
    {
        return fragmentEdges - fragmentEdges / duplicateSize -
            (fragmentEdges % duplicateSize != 0);
    }

    /**
     * @brief Old branch and bound heuristic, used only during initial enumeration
     * 
     * @return int (maximum duplicatable bonds value). Lower bound MA is total
     * bonds - 1 - this value.
     */
    int maxDupBonds()
    {
        return maxDupBonds(edgeCounts[0]);
    }

    /**
     * @brief Bound duplicatable bonds using eligible edges in each fragment
     * 
     * @param targetMasks The vector of bitsets used in place of sizeList
     */
    int maxDupBonds(
        int maxFragSize,
        const vector<EdgeMask> &targetMasks
    ) const
    {
        int dupBondsTotal = -ceilLog2(maxFragSize);
        for (size_t i = 0; i < masks.size(); i++)
        {
            dupBondsTotal += fixedSizeDupBondsForFragment(
                edgeCounts[i],
                maxFragSize,
                static_cast<int>(targetMasks[i].count())
            );
        }
        return dupBondsTotal;
    }

    /**
     * @brief Like the function above but finds the maximum duplicate bonds for a vector of vector of bitsets
     * 
     * @param fragSizeListMax The result vector
     * @param targetMasks The vector of vector of bitsets
     */
    void maxDupBondsPrefix(
        vi &fragSizeListMax,
        int maxFragSize,
        const vector<vector<EdgeMask> > &targetMasks
    ) const
    {
        if (maxFragSize < 2)
        {
            fragSizeListMax.clear();
            return;
        }

        fragSizeListMax.assign(maxFragSize - 1, 0);
        for (const int edgeCount : edgeCounts)
        {
            fragSizeListMax[0] += edgeCount / 2;
        }
        fragSizeListMax[0]--;

        for (int duplicateSize = 3;
             duplicateSize <= maxFragSize;
             duplicateSize++)
        {
            const size_t index = duplicateSize - 2;
            int bound = -ceilLog2(duplicateSize);
            const vector<EdgeMask> &duplicateMasks = targetMasks[index];
            for (size_t i = 0; i < masks.size(); i++)
            {
                bound += fixedSizeDupBondsForFragment(
                    edgeCounts[i],
                    duplicateSize,
                    static_cast<int>(duplicateMasks[i].count())
                );
            }
            fragSizeListMax[index] = max(
                bound,
                fragSizeListMax[index - 1]
            );
        }
    }

    /**
     * @brief The simple branch and bound from v4
     * 
     * @param maxFragSize The maximum allowed fragment size
     * @return int The maximum number of duplicatable bonds
     */
    int maxDupBonds(int maxFragSize) const
    {
        int dupBonds2 = 0, dupBondsTotal;

        for (const int edgeCount : edgeCounts) dupBonds2 += edgeCount / 2;
        dupBonds2--;
        for (int duplicateSize = 3;
             duplicateSize <= maxFragSize;
             duplicateSize++)
        {
            dupBondsTotal = 0;
            for (const int edgeCount : edgeCounts)
            {
                dupBondsTotal += unrestrictedDupBondsForFragment(
                    edgeCount,
                    duplicateSize
                );
            }
            dupBondsTotal -= ceilLog2(duplicateSize);
            if (dupBondsTotal > dupBonds2) dupBonds2 = dupBondsTotal;
        }
        return dupBonds2;
    }

    /**
     * @brief Old branch and bound. Only used during initial enumeration
     * 
     * @return int The Lower bound
     */
    int lowBoundAI()
    {
        return static_cast<int>(totalBonds) - sumDupBonds - 1 - maxDupBonds();
    }

    /**
     * @brief Upper bound MA given the sum of duplicatable bonds
     * 
     * @return int The upper bound
     */
    int AI()
    {
        return static_cast<int>(totalBonds) - sumDupBonds - 1;
    }

    /**
     * @brief Build the canonical fragment key used by the transposition table.
     *
     * @param sorted Reused storage populated with the canonical key
     */
    void assemblyHashCalculator(vi &sorted)
    {
        sorted.assign(masks.size(), -1);
        for (size_t i = 0; i < masks.size(); i++)
        {
            const auto entry = bitsetHashTable.find(masks[i]);
            if (entry != bitsetHashTable.end()) sorted[i] = entry->second.first;
        }
        sort(sorted.begin() + 1, sorted.end());
    }
};
