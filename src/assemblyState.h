#pragma once

#include "additionChainBounds.h"

/**
 * @brief Assembly state data structure. Records the current state of this assembly pathway
 */
struct assemblyState
{
    /// @brief Fragment masks and metadata kept valid as one unit
    vector<assemblyFragment> fragments;
    /// @brief Join savings accumulated by removing duplicate occurrences
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

    void appendFragment(
        EdgeMask &&mask,
        int edgeCount,
        int canonicalId = unknownCanonicalId,
        bool connected = false
    )
    {
        fragments.emplace_back(std::move(mask), edgeCount, canonicalId, connected);
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
     * @brief Bound further join savings using the first fragment's size limit
     *
     * @pre The state contains a first fragment defining the duplicate-size limit.
     * @return Upper bound on additional savings. Subtract this value from
     * assemblyIndex() to obtain a lower bound on the assembly index.
     */
    int maxDupBonds() const
    {
        return maxDupBonds(fragments[0].edgeCount);
    }

    /**
     * @brief Bound duplicatable bonds using eligible edges in each fragment
     *
     * @param maximumFragmentSize Duplicate size, at least two
     * @param targetMasks Eligible-edge mask or count for each state fragment
     */
    template<typename MaskRange>
    int maxDupBonds(
        int maximumFragmentSize,
        const MaskRange &targetMasks
    ) const
    {
        int totalDuplicateBondBound =
            -assembly_bounds::scalarLowerBound(maximumFragmentSize);
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
     * @brief Build prefix bounds on further savings for increasing duplicate sizes
     *
     * @param maximumByFragmentSize Output indexed by duplicate size minus two;
     * each entry bounds savings for that size or any smaller size
     * @param maximumFragmentSize Largest duplicate size to consider
     * @param targetMasks Eligible-edge counts or masks indexed by duplicate
     * size minus two, then by state fragment; size two uses the unrestricted bound
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
            int bound = -assembly_bounds::scalarLowerBound(duplicateSize);
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
     * @brief Bound further join savings without restricting eligible edges
     *
     * @param maximumFragmentSize The maximum allowed fragment size
     * @return Upper bound on additional savings across permitted duplicate sizes
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
            candidateDuplicateBondBound -=
                assembly_bounds::scalarLowerBound(duplicateSize);
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
     * @brief Bound the assembly index using unrestricted duplicate eligibility
     *
     * @return Lower bound on the assembly index
     */
    int lowerBoundAssemblyIndex() const
    {
        return max(
            assemblyCompositionLowerBound,
            static_cast<int>(totalBonds) - sumDupBonds - 1 - maxDupBonds()
        );
    }

    /**
     * @brief Upper bound on the assembly index from savings already realised
     *
     * @return int The upper bound
     */
    int assemblyIndex() const
    {
        return static_cast<int>(totalBonds) - sumDupBonds - 1;
    }

    /**
     * @brief Build the canonical fragment key used by the transposition table.
     *
     * The retained fragment stays first; the remaining IDs are sorted. Every
     * fragment must already have a resolved canonical ID.
     *
     * @param sorted Reused storage populated with the canonical key
     */
    void assemblyHashCalculator(IntegerVector &sorted) const
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
