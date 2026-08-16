/**
 * @brief vector<int> hash, used to hash assembly states
 */
template<>
struct std::hash<vi>
{
    size_t operator()(const vi &v) const
    {
        std::size_t seed = v.size();
        for(auto& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

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
    assemblyPath * ap;

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
    /// @brief index of the state
    int ix = 0;
    /// @brief path that was used to generate this state
    assemblyPath * apPtr = nullptr;

   /**
   * @brief Return the maximum fragment size, by counting the number of set bits in the first mask
   *
   * @return int (maximum size of a fragment)
   */
    int maxFragSizeF()
    {
        return masks[0].count();
    }

    /**
     * @brief Old branch and bound heuristic, used only during initial enumeration
     * 
     * @return int (maximum duplicatable bonds value). Lower bound MA is total
     * bonds - 1 - this value.
     */
    int maxDupBonds()
    {
        int dupBonds2 = 0, dupBondsTotal, maxFragSize = maxFragSizeF();
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
            dupBondsTotal -= ceil(log2(j));
            if (dupBondsTotal > dupBonds2) dupBonds2 = dupBondsTotal;
        }
        return dupBonds2;
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
        int dupBonds2 = 0, dupBondsTotal;

        int j = maxFragSize;
        vi adjustedSizeList(sizeList.size()), adjustedSizeList2(sizeList.size());
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            adjustedSizeList[i] = sizeList[i] - sizeList[i] % j;
            adjustedSizeList2[i] = sizeListMain[i] - adjustedSizeList[i];
        }
        dupBondsTotal = 0;
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            dupBondsTotal += (adjustedSizeList[i] - adjustedSizeList[i]/j);
            dupBondsTotal += (adjustedSizeList2[i] - adjustedSizeList2[i]/(j - 1));
            if (adjustedSizeList2[i] % (j - 1) != 0) dupBondsTotal--;
        }
        dupBondsTotal -= ceil(log2(j));
        return dupBondsTotal;
    }

    /**
     * @brief Basically same function as above but takes a vector of bitsets and uses the count to find sizeList
     * 
     * @param targetMasks The vector of bitsets used in place of sizeList
     */
    int maxDupBonds(vi &sizeListMain, int maxFragSize, vector<standardBitset> &targetMasks)
    {
        int dupBonds2 = 0, dupBondsTotal;
        vi sizeList(targetMasks.size());
        
        for (size_t i = 0; i < masks.size(); i++)
        {
            sizeList[i] = targetMasks[i].count();
        }

        int j = maxFragSize;
        vi adjustedSizeList(sizeList.size()), adjustedSizeList2(sizeList.size());
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            adjustedSizeList[i] = sizeList[i] - sizeList[i] % j;
            adjustedSizeList2[i] = sizeListMain[i] - adjustedSizeList[i];
        }
        dupBondsTotal = 0;
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            dupBondsTotal += (adjustedSizeList[i] - adjustedSizeList[i]/j);
            dupBondsTotal += (adjustedSizeList2[i] - adjustedSizeList2[i]/(j - 1));
            if (adjustedSizeList2[i] % (j - 1) != 0) dupBondsTotal--;
        }
        dupBondsTotal -= ceil(log2(j));
        return dupBondsTotal;
    }

    /**
     * @brief Like the function above but finds the maximum duplicate bonds for a vector of vector of bitsets
     * 
     * @param fragSizeList The result vector
     * @param targetMasks The vector of vector of bitsets
     */
    void maxDupBonds(vi &fragSizeList, int maxFragSize, vector<vector<standardBitset> > &targetMasks)
    {
        int dupBonds2 = 0, dupBondsTotal;
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
            vi &sizeList = sizeLists[j - 2], &sizeList2 = sizeLists[0];
            vi adjustedSizeList(sizeList.size()), adjustedSizeList2(sizeList.size());
            for (size_t i = 0; i < sizeList.size(); i++)
            {
                adjustedSizeList[i] = sizeList[i] - sizeList[i] % j;
                adjustedSizeList2[i] = sizeList2[i] - adjustedSizeList[i];
            }
            dupBondsTotal = 0;
            for (size_t i = 0; i < sizeList.size(); i++)
            {
                dupBondsTotal += (adjustedSizeList[i] - adjustedSizeList[i]/j);
                dupBondsTotal += (adjustedSizeList2[i] - adjustedSizeList2[i]/(j - 1));
                if (adjustedSizeList2[i] % (j - 1) != 0) dupBondsTotal--;
            }
            dupBondsTotal -= ceil(log2(j));
            fragSizeList[j - 2] = dupBondsTotal;
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
            dupBondsTotal -= ceil(log2(j));
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
     * @brief Calculates the hash for the assembly state by returning a vi which is subjected to the vi hash
     * 
     * @return vi The vector<int> to be hashed
     */
    vi assemblyHashCalculator()
    {
        vi sorted(masks.size(), -1);
        for (size_t i = 0; i < masks.size(); i++)
        {
            if (bitsetHashTable.count(masks[i])) sorted[i] = bitsetHashTable[masks[i]].first;
        }
        sort(sorted.begin() + 1, sorted.end());
        return sorted;
    }
    
    /**
     * @brief Prints the assembly state
     * 
     */
    void print()
    {
        cout << "printing assembly state: " << ix << '\n';
        cout << "masks:\n";
        for (size_t i = 0; i < masks.size(); i++)
        {
            cout << masks[i] << '\n';
        }
        for (size_t i = 0; i < masks.size(); i++)
        {
            cout << "Fragment: " << i <<'\n';
            bool isCyclic;
            molGraph mg = constructFromEdgeList(targetMolecule, univEdgeList, masks[i], isCyclic);
            mg.printToCout();
        }
        cout << "AI: " << AI() << '\n';
        cout << "lowbound: " << lowBoundAI() << '\n';
    }
};