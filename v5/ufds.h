/**
 * @brief node of a disjoint set
 */
struct disjointSetNode
{
    int parent = -1, rank = 0;
    disjointSetNode(){}
    disjointSetNode(int _parent): parent(_parent){}
};

/**
 * @brief Disjoint set data structure for constructing an edge list from a bitmask.
 * Practically identical to textbook UFDS data structure
 */
struct disjointSet
{
    /// nodes
    vector<disjointSetNode> elements;
    
    disjointSet(size_t size){elements.resize(size);}

    /// standard disjoint set function
    size_t find(size_t idx)
    {
        if (elements[idx].parent != idx)
        {
            elements[idx].parent = find(elements[idx].parent);
        }
        return elements[idx].parent;
    }

    /// standard disjoint set function
    void insert(int target, int parent)
    {
        elements[target].parent = parent;
    }

    /// standard disjoint set function
    bool merge(size_t x, size_t y)
    {
        size_t rootx = find(x), rooty = find(y);
        if (rootx == rooty) return true;
        if (elements[rootx].rank > elements[rooty].rank)
        {
            elements[rooty].parent = rootx;
        }
        else
        {
            elements[rootx].parent = rooty;
            if (elements[rootx].rank == elements[rooty].rank) elements[rooty].rank++;
        }
        return false;
    }
};
 
/**
 * @brief for UFDS split node - variant on textbook UFDS
 */
struct ufdsSplitNode
{
    /// parent, rank, fragment this is part of
    int parent = -1, rank = 0, val = -1;

    ufdsSplitNode(){}
    ufdsSplitNode(int _parent, int _val): parent(_parent)
    {
        val = _val;
    }
};

/**
 * @brief for UFDS split node - variant on textbook UFDS
 */
struct ufdsSplit
{
    vector<ufdsSplitNode> elements;
    vector<pii> extraVals;
    int maxElement = 0;
    
    ufdsSplit(size_t size){elements.resize(size);}

    /// standard disjoint set function
    size_t find(size_t idx)
    {
        if (elements[idx].parent != idx)
        {
            elements[idx].parent = find(elements[idx].parent);
        }
        return elements[idx].parent;
    }

    /**
     * @brief Used if one atom has not been seen before
     *
     */
    void insert(int target, int parent, int val)
    {
        if (target > maxElement) maxElement = target;
        if (parent > maxElement) maxElement = parent;
        ufdsSplitNode u(parent, val);
        elements[target] = u;
    }
    
    /**
     * @brief Used if both atoms have not been seen before
     *
     */
    void doubleInsert(int target, int parent, int val)
    {
        if (target > maxElement) maxElement = target;
        if (parent > maxElement) maxElement = parent;
        ufdsSplitNode u(parent, val);
        elements[target] = u;
        elements[parent] = u;
    }

    /**
     * @brief Used if both atoms have been seen before
     *
     */
    void merge(size_t x, size_t y, int yval)
    {
        size_t rootx = find(x), rooty = find(y);
        extraVals.push_back(pii(rooty, yval));
        if (rootx != rooty)
        {
            if (elements[rootx].rank > elements[rooty].rank)
            {
                elements[rooty].parent = rootx;
            }
            else
            {
                elements[rootx].parent = rooty;
                if (elements[rootx].rank == elements[rooty].rank) elements[rooty].rank++;
            }
        }
    }

    /**
     * @brief The splitting function used during the fragmentation
     *
     * @param maskList The output
     */
    void split(vector<standardBitset> &maskList)
    {
        vi uniques(maxElement + 1, -1);
        vector<standardBitset> tempMaskList;
        for (size_t i = 0; i <= maxElement; i++)
        {
            if (elements[i].parent != -1)
            {
                find(i);
                if (uniques[elements[i].parent] == -1)
                {
                    uniques[elements[i].parent] = tempMaskList.size();
                    standardBitset b = 0; b.set(elements[i].val);
                    tempMaskList.push_back(b);
                }
                else
                {
                    tempMaskList[uniques[elements[i].parent]][elements[i].val] = 1;
                }
            }
        }
        for (size_t i = 0; i < extraVals.size(); i++)
        {
            tempMaskList[uniques[elements[extraVals[i].first].parent]][extraVals[i].second] = 1;
        }
        for (size_t i = 0; i < tempMaskList.size(); i++)
        {
            if (tempMaskList[i].count() > 1)
            {
                maskList.push_back(tempMaskList[i]);
            }
        }
    }
};


