/**
 * @brief Return the low machine word of a bitset, masked below a limit
 */
template<typename Bitset>
unsigned long long bitsetLowWordBelow(const Bitset &mask, size_t limit)
{
    constexpr size_t wordBits = numeric_limits<unsigned long long>::digits;
    const Bitset lowWordMask(numeric_limits<unsigned long long>::max());
    unsigned long long word = (mask & lowWordMask).to_ullong();
    if (limit < wordBits) word &= (1ULL << limit) - 1;
    return word;
}

/**
 * @brief Visit set bit indices below a limit in ascending order
 *
 * Sparse implementations visit only set bits; dense masks and the portable
 * fallback use a bounded linear scan.
 */
template<typename Bitset, typename Visitor>
void forEachSetBitWithWideLimit(
    const Bitset &mask,
    size_t limit,
    Visitor &&visitor
)
{
    limit = min(limit, mask.size());
    if constexpr (requires(const Bitset &bits, size_t index) {
        bits._Find_first();
        bits._Find_next(index);
    })
    {
        // A linear scan is cheaper for dense residuals and is still
        // proportional to them when at least one bit in four is set.
        if (mask.count() * 4 >= limit)
        {
            for (size_t index = 0; index < limit; index++)
            {
                if (mask[index]) visitor(index);
            }
            return;
        }
        for (
            size_t index = mask._Find_first();
            index < limit;
            index = mask._Find_next(index)
        )
        {
            visitor(index);
        }
    }
    else
    {
        for (size_t index = 0; index < limit; index++)
        {
            if (mask[index]) visitor(index);
        }
    }
}

template<typename Bitset, typename Visitor>
void forEachSetBitBelow(const Bitset &mask, size_t limit, Visitor &&visitor)
{
    limit = min(limit, mask.size());
    constexpr size_t wordBits = numeric_limits<unsigned long long>::digits;
    if (limit > wordBits)
    {
        forEachSetBitWithWideLimit(
            mask,
            limit,
            std::forward<Visitor>(visitor)
        );
        return;
    }
    unsigned long long word = bitsetLowWordBelow(mask, limit);
    while (word != 0)
    {
        const size_t index = std::countr_zero(word);
        visitor(index);
        word &= word - 1;
    }
}

/**
 * @brief node of a disjoint set
 */
struct disjointSetNode
{
    int parent = -1, rank = 0;
    disjointSetNode() = default;
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
        if (elements[idx].parent != static_cast<int>(idx))
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
    int16_t parent = -1, rank = 0, val = -1, component = -1;
    uint32_t generation = 0;

    ufdsSplitNode() = default;
    ufdsSplitNode(int _parent, int _val, uint32_t _generation):
        parent(static_cast<int16_t>(_parent)),
        val(static_cast<int16_t>(_val)),
        generation(_generation)
    {
    }
};

static_assert(BITSET_LENGTH <= numeric_limits<int16_t>::max());

/**
 * @brief for UFDS split node - variant on textbook UFDS
 */
struct ufdsSplit
{
    static constexpr size_t atomWordBits = numeric_limits<uint64_t>::digits;
    static constexpr size_t atomWordCount =
        (BITSET_LENGTH + atomWordBits - 1) / atomWordBits;
    static_assert(atomWordCount <= numeric_limits<uint64_t>::digits);

    vector<ufdsSplitNode> elements;
    vector<pii> extraVals;
    array<uint64_t, atomWordCount> touchedAtomWords{};
    uint64_t touchedAtomWordMask = 0;
    // Keep the active generation distinct from default-constructed nodes even
    // before the first reset.
    uint32_t generation = 1;
    
    void reset()
    {
        generation++;
        if (generation == 0)
        {
            for (ufdsSplitNode &element : elements) element.generation = 0;
            generation = 1;
        }
        if ((touchedAtomWordMask & 1) != 0) touchedAtomWords[0] = 0;
        uint64_t activeWords = touchedAtomWordMask & ~uint64_t{1};
        while (activeWords != 0)
        {
            const size_t wordIndex = std::countr_zero(activeWords);
            touchedAtomWords[wordIndex] = 0;
            activeWords &= activeWords - 1;
        }
        touchedAtomWordMask = 0;
        extraVals.clear();
    }

    void markTouched(size_t index)
    {
        const size_t wordIndex = index / atomWordBits;
        touchedAtomWords[wordIndex] |= uint64_t{1} << (index % atomWordBits);
        touchedAtomWordMask |= uint64_t{1} << wordIndex;
    }

    bool contains(size_t index) const
    {
        return elements[index].generation == generation;
    }

    /// standard disjoint set function
    size_t find(size_t idx)
    {
        if (elements[idx].parent != static_cast<int>(idx))
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
        markTouched(target);
        elements[target] = ufdsSplitNode(parent, val, generation);
    }
    
    /**
     * @brief Used if both atoms have not been seen before
     *
     */
    void doubleInsert(int target, int parent, int val)
    {
        markTouched(target);
        if (parent != target)
        {
            markTouched(parent);
        }
        ufdsSplitNode u(parent, val, generation);
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
     * @param tempMaskList Reusable component-mask buffer; must not alias maskList
     */
    void splitWithBuffers(
        vector<standardBitset> &maskList,
        vector<standardBitset> &tempMaskList
    )
    {
        tempMaskList.clear();
        auto addTouchedAtom = [&](size_t index) {
            find(index);
            const int root = elements[index].parent;
            int16_t &component = elements[root].component;
            if (component == -1)
            {
                component = tempMaskList.size();
                standardBitset b = 0; b.set(elements[index].val);
                tempMaskList.push_back(b);
            }
            else
            {
                tempMaskList[component].set(elements[index].val);
            }
        };
        uint64_t atoms = touchedAtomWords[0];
        while (atoms != 0)
        {
            const size_t bitIndex = std::countr_zero(atoms);
            addTouchedAtom(bitIndex);
            atoms &= atoms - 1;
        }
        uint64_t activeWords = touchedAtomWordMask & ~uint64_t{1};
        while (activeWords != 0)
        {
            const size_t wordIndex = std::countr_zero(activeWords);
            atoms = touchedAtomWords[wordIndex];
            while (atoms != 0)
            {
                const size_t bitIndex = std::countr_zero(atoms);
                addTouchedAtom(wordIndex * atomWordBits + bitIndex);
                atoms &= atoms - 1;
            }
            activeWords &= activeWords - 1;
        }
        for (size_t i = 0; i < extraVals.size(); i++)
        {
            const size_t root = find(extraVals[i].first);
            const int component = elements[root].component;
            tempMaskList[component].set(extraVals[i].second);
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
