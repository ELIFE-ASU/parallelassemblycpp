#pragma once

/**
 * @brief Assembly state data structure. Records the current state of this assembly pathway
 */
struct assemblyState
{
    /// @brief Fragment masks and metadata kept valid as one unit
    vector<assemblyFragment> fragments;
    /// @brief number of duplicated bonds
    int sumDupBonds = 0;

    void appendFragment(
        const EdgeMask &mask,
        int edgeCount,
        int canonicalId = unknownCanonicalId,
        bool connected = false
    )
    {
        fragments.emplace_back(mask, edgeCount, canonicalId, connected);
    }

    void appendFragment(const assemblyFragment &fragment)
    {
        fragments.push_back(fragment);
    }

    void reserveFragments(size_t capacity)
    {
        fragments.reserve(capacity);
    }

    void clearFragments()
    {
        fragments.clear();
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
        return maxDupBonds(fragments[0].edgeCount);
    }

    /**
     * @brief Bound duplicatable bonds using eligible edges in each fragment
     *
     * @param targetMasks The vector of bitsets used in place of sizeList
     */
    template<typename MaskRange>
    int maxDupBonds(
        int maximumFragmentSize,
        const MaskRange &targetMasks
    ) const
    {
        int totalDuplicateBondBound = -ceilLog2(maximumFragmentSize);
        for (
            size_t fragmentIndex = 0;
            fragmentIndex < fragments.size();
            fragmentIndex++
        )
        {
            totalDuplicateBondBound += fixedSizeDupBondsForFragment(
                fragments[fragmentIndex].edgeCount,
                maximumFragmentSize,
                static_cast<int>(maskCountAt(targetMasks, fragmentIndex))
            );
        }
        return totalDuplicateBondBound;
    }

    /**
     * @brief Like the function above but finds the maximum duplicate bonds for a vector of vector of bitsets
     *
     * @param maximumByFragmentSize The result vector
     * @param targetMasks The vector of vector of bitsets
     */
    template<typename MaskTable>
    void maxDupBondsPrefix(
        IntegerVector &maximumByFragmentSize,
        int maximumFragmentSize,
        const MaskTable &targetMasks
    ) const
    {
        if (maximumFragmentSize < 2)
        {
            maximumByFragmentSize.clear();
            return;
        }

        maximumByFragmentSize.assign(maximumFragmentSize - 1, 0);
        for (const assemblyFragment &fragment : fragments)
        {
            maximumByFragmentSize[0] += fragment.edgeCount / 2;
        }
        maximumByFragmentSize[0]--;

        for (int duplicateSize = 3;
             duplicateSize <= maximumFragmentSize;
             duplicateSize++)
        {
            const size_t index = duplicateSize - 2;
            int bound = -ceilLog2(duplicateSize);
            for (
                size_t fragmentIndex = 0;
                fragmentIndex < fragments.size();
                fragmentIndex++
            )
            {
                bound += fixedSizeDupBondsForFragment(
                    fragments[fragmentIndex].edgeCount,
                    duplicateSize,
                    static_cast<int>(tableMaskCountAt(
                        targetMasks,
                        index,
                        fragmentIndex
                    ))
                );
            }
            maximumByFragmentSize[index] = max(
                bound,
                maximumByFragmentSize[index - 1]
            );
        }
    }

    /**
     * @brief The simple branch and bound from v4
     *
     * @param maximumFragmentSize The maximum allowed fragment size
     * @return int The maximum number of duplicatable bonds
     */
    int maxDupBonds(int maximumFragmentSize) const
    {
        int bestDuplicateBondBound = 0;

        for (const assemblyFragment &fragment : fragments)
            bestDuplicateBondBound += fragment.edgeCount / 2;
        bestDuplicateBondBound--;
        for (int duplicateSize = 3;
             duplicateSize <= maximumFragmentSize;
             duplicateSize++)
        {
            int candidateDuplicateBondBound = 0;
            for (const assemblyFragment &fragment : fragments)
            {
                candidateDuplicateBondBound += unrestrictedDupBondsForFragment(
                    fragment.edgeCount,
                    duplicateSize
                );
            }
            candidateDuplicateBondBound -= ceilLog2(duplicateSize);
            bestDuplicateBondBound = max(
                bestDuplicateBondBound,
                candidateDuplicateBondBound
            );
        }
        return bestDuplicateBondBound;
    }

private:
    template<typename MaskRange>
    static size_t maskCountAt(const MaskRange &masks, size_t index)
    {
        if constexpr (requires { masks.maskCount(index); })
            return masks.maskCount(index);
        else
            return masks[index].count();
    }

    template<typename MaskTable>
    static size_t tableMaskCountAt(
        const MaskTable &masks,
        size_t row,
        size_t column
    )
    {
        if constexpr (requires { masks.maskCount(row, column); })
            return masks.maskCount(row, column);
        else
            return masks[row][column].count();
    }

public:

    /**
     * @brief Old branch and bound. Only used during initial enumeration
     *
     * @return int The Lower bound
     */
    int lowerBoundAssemblyIndex()
    {
        return static_cast<int>(totalBonds) - sumDupBonds - 1 - maxDupBonds();
    }

    /**
     * @brief Upper bound MA given the sum of duplicatable bonds
     *
     * @return int The upper bound
     */
    int assemblyIndex()
    {
        return static_cast<int>(totalBonds) - sumDupBonds - 1;
    }

    /**
     * @brief Build the canonical fragment key used by the transposition table.
     *
     * @param sorted Reused storage populated with the canonical key
     */
    void assemblyHashCalculator(IntegerVector &sorted)
    {
        sorted.resize(fragments.size());
        for (
            size_t fragmentIndex = 0;
            fragmentIndex < fragments.size();
            fragmentIndex++
        )
            sorted[fragmentIndex] = fragments[fragmentIndex].canonicalId;
        if (sorted.size() > 1) sort(sorted.begin() + 1, sorted.end());
    }
};
