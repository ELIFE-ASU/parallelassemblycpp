/**
 * @brief Struct which is inserted into the hash table to store assembly states and recover the pathway
 */
struct assemblyPath
{
    /// vi represneting canonical indices of the fragments of the assembly state
    vi key;
    /// @brief number of duplicated bonds found so far
    int sumDupBonds;
    /// @brief needed to reconstruct the pathway
    unsigned short match, duplicate;
    /// @brief Assembly state from which this state is generated. needed to reconstruct the pathway
    assemblyPath * parent;
};

/// Pointer for the minimum assembly path
assemblyPath * minAssemblyPath = nullptr;

/**
 * @brief Wrapper for pathway hash table because C++ unordered_map does not guarantee pointer will remain unchanged
 * 
 */
struct apWrapper
{
    mutable assemblyPath * ap;

    bool operator == (const apWrapper&ap2) const
    {
        return ap->key == ap2.ap->key;
    }
};

/**
 * @brief vi hash called by apWrapper
 * 
 */
template<>
struct std::hash<apWrapper>
{
    size_t operator()(const apWrapper &ap) const
    {
        std::size_t seed = ap.ap->key.size();
        for(auto& i : ap.ap->key) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

/// Hash table for assembly states for pathway algorithm
std::unordered_set<apWrapper> pathAssemblyMap;

/**
 * @brief Assembly state data structure. Records the current state of this assembly pathway
 */
struct assemblyState
{   
    /// @brief each mask represents a separate fragment as a boolean edge list
    vector<standardBitset> masks;
    /// @brief number of duplicated bonds
    int sumDupBonds = 0;
    /// @brief path that was used to generate this state
    assemblyPath * apPtr = nullptr;

    /**
     * @brief Old branch and bound heuristic, used only during initial enumeration
     * 
     * @return int (maximum duplicatable bonds value). Lower bound MA is total
     * bonds - 1 - this value.
     */
    int maxDupBonds()
    {
        return maxDupBonds(static_cast<int>(masks[0].count()));
    }

    /**
     * @brief Calculates the maximum number of duplicatable bonds if given a maximum fragment size and a maximal fragment list
     * 
     * @param sizeListMain The bitset counts of each fragment
     * @param maxFragSize The maximum allowed fragment size
     * @param sizeList The bitset counts of all duplicates in each fragment
     * @return int The maximum number of duplicatable bonds
     */
    int maxDupBonds(vi &sizeListMain, int maxFragSize, vi &sizeList)
    {
        const int j = maxFragSize;
        int dupBondsTotal = 0;
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            const int adjustedSize = sizeList[i] - sizeList[i] % j;
            const int remainingSize = sizeListMain[i] - adjustedSize;
            dupBondsTotal += adjustedSize - adjustedSize / j;
            dupBondsTotal += remainingSize - remainingSize / (j - 1);
            if (remainingSize % (j - 1) != 0) dupBondsTotal--;
        }
        dupBondsTotal -= ceilLog2(j);
        return dupBondsTotal;
    }

    /**
     * @brief Basically same function as above but takes a vector of bitsets and uses the count to find sizeList
     * 
     * @param targetMasks The vector of bitsets used in place of sizeList
     */
    int maxDupBonds(vi &sizeListMain, int maxFragSize, vector<standardBitset> &targetMasks)
    {
        vi sizeList(targetMasks.size());
        
        for (size_t i = 0; i < masks.size(); i++)
        {
            sizeList[i] = static_cast<int>(targetMasks[i].count());
        }
        return maxDupBonds(sizeListMain, maxFragSize, sizeList);
    }

    /**
     * @brief Like the function above but finds the maximum duplicate bonds for a vector of vector of bitsets
     * 
     * @param fragSizeList The result vector
     * @param targetMasks The vector of vector of bitsets
     */
    void maxDupBonds(vi &fragSizeList, int maxFragSize, vector<vector<standardBitset> > &targetMasks)
    {
        int dupBonds2 = 0;
        fragSizeList.resize(maxFragSize - 1);
        vector<vi> sizeLists(fragSizeList.size());
        
        for (size_t j = 0; j < targetMasks.size(); j++)
        {
            sizeLists[j].assign(targetMasks[j].size(), 0);
            for (size_t i = 0; i < masks.size(); i++)
            {
                sizeLists[j][i] = targetMasks[j][i].count();
            }
        }

        for (size_t i = 0; i < sizeLists[0].size(); i++) dupBonds2 += sizeLists[0][i]/2;
        dupBonds2--;
        fragSizeList[0] = dupBonds2;

        for (int j = 3; j <= maxFragSize; j++)
        {
            fragSizeList[j - 2] = maxDupBonds(
                sizeLists[0],
                j,
                sizeLists[j - 2]
            );
        }
    }

    /**
     * @brief The simple branch and bound from v4
     * 
     * @param maxFragSize The maximum allowed fragment size
     * @return int The maximum number of duplicatable bonds
     */
    int maxDupBonds(int maxFragSize)
    {
        int dupBonds2 = 0, dupBondsTotal;
        vi sizeList(masks.size());
        
        for (size_t i = 0; i < masks.size(); i++)
        {
            sizeList[i] = masks[i].count();
        }

        for (size_t i = 0; i < sizeList.size(); i++) dupBonds2 += sizeList[i]/2;
        dupBonds2--;
        for (int j = 3; j <= maxFragSize; j++)
        {
            dupBondsTotal = 0;
            for (size_t i = 0; i < sizeList.size(); i++)
            {
                dupBondsTotal += (sizeList[i] - sizeList[i]/j);
                if (sizeList[i] % j != 0) dupBondsTotal--;
            }
            dupBondsTotal -= ceilLog2(j);
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
        return totalBonds - sumDupBonds - 1 - maxDupBonds();
    }

    /**
     * @brief Function for calculating lower bound on MA given an estimate and a maximum allowed fragment size
     * 
     * @param maxFragSize The maximum allowed fragment size
     * @param estimate The estimate
     * @return int The lower bound
     */
    int lowBoundAI(int maxFragSize, int estimate)
    {
        if (maxFragSize > 2) estimate = max(estimate, maxDupBonds(maxFragSize - 1));
        return totalBonds - sumDupBonds - 1 - estimate;
    }

    /**
     * @brief Upper bound MA given the sum of duplicatable bonds
     * 
     * @return int The upper bound
     */
    int AI()
    {
        return totalBonds - sumDupBonds - 1;
    }

    /**
     * @brief Build the canonical fragment key used by apWrapper's hash
     *
     * @return vi The vector<int> to be hashed
     */
    vi assemblyHashCalculator()
    {
        vi sorted(masks.size(), -1);
        for (size_t i = 0; i < masks.size(); i++)
        {
            const auto entry = bitsetHashTable.find(masks[i]);
            if (entry != bitsetHashTable.end()) sorted[i] = entry->second.first;
        }
        sort(sorted.begin() + 1, sorted.end());
        return sorted;
    }
};
