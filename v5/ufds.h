/**
 * @brief Return the low machine word of a bitset, masked below a limit
 */
template<typename Bitset>
unsigned long long bitsetLowWordBelow(const Bitset &mask, size_t limit)
{
    if constexpr (requires { mask.lowWordBelow(limit); })
    {
        return mask.lowWordBelow(limit);
    }
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
    if constexpr (requires(const Bitset &bits, size_t wordIndex) {
        Bitset::activeWordCount();
        bits.activeWord(wordIndex);
    })
    {
        constexpr size_t wordBits =
            numeric_limits<unsigned long long>::digits;
        const size_t wordCount =
            limit / wordBits + (limit % wordBits != 0);
        for (size_t wordIndex = 0; wordIndex < wordCount; wordIndex++)
        {
            unsigned long long word = mask.activeWord(wordIndex);
            if (wordIndex + 1 == wordCount && limit % wordBits != 0)
            {
                word &= (1ULL << (limit % wordBits)) - 1;
            }
            while (word != 0)
            {
                visitor(
                    wordIndex * wordBits +
                    static_cast<size_t>(std::countr_zero(word))
                );
                word &= word - 1;
            }
        }
        return;
    }
    if constexpr (requires(const Bitset &bits, size_t index) {
        bits.findFirst();
        bits.findNext(index);
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
            size_t index = mask.findFirst();
            index < limit;
            index = mask.findNext(index)
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
    int16_t parent = -1, rank = 0;
    int32_t val = -1, component = -1;
    uint32_t generation = 0;

    ufdsSplitNode() = default;
    ufdsSplitNode(int _parent, int _val, uint32_t _generation):
        parent(static_cast<int16_t>(_parent)),
        val(_val),
        generation(_generation)
    {
    }
};

static_assert(sizeof(ufdsSplitNode) == 16);

/**
 * @brief for UFDS split node - variant on textbook UFDS
 */
struct ufdsSplit
{
    static constexpr size_t atomWordBits = numeric_limits<uint64_t>::digits;
    // Keep the common case allocation-free and identical to the former
    // fixed-capacity representation. Larger atom indices spill into the
    // vectors below, removing this structure's former 512-atom boundary;
    // molGraph's adjacency-index type retains a separate, much higher limit.
    static constexpr size_t inlineAtomBitCount = 512;
    static constexpr size_t inlineAtomWordCount =
        inlineAtomBitCount / atomWordBits;

    vector<ufdsSplitNode> elements;
    vector<pii> extraVals;
    array<uint64_t, inlineAtomWordCount> touchedAtomWords{};
    uint64_t touchedAtomWordMask = 0;
    vector<uint64_t> wideTouchedAtomWords;
    vector<size_t> touchedWideAtomWordIndices;
    vector<uint64_t> componentMaskWords;
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
        if (!touchedWideAtomWordIndices.empty()) [[unlikely]]
        {
            for (const size_t wordIndex : touchedWideAtomWordIndices)
            {
                wideTouchedAtomWords[wordIndex - inlineAtomWordCount] = 0;
            }
            touchedWideAtomWordIndices.clear();
        }
        extraVals.clear();
    }

    [[gnu::always_inline]] void markTouched(size_t index)
    {
        const size_t wordIndex = index / atomWordBits;
        const uint64_t bit = uint64_t{1} << (index % atomWordBits);
        if (wordIndex < inlineAtomWordCount) [[likely]]
        {
            touchedAtomWords[wordIndex] |= bit;
            touchedAtomWordMask |= uint64_t{1} << wordIndex;
            return;
        }

        markTouchedWide(wordIndex, bit);
    }

    [[gnu::noinline]] void markTouchedWide(size_t wordIndex, uint64_t bit)
    {
        const size_t wideWordIndex = wordIndex - inlineAtomWordCount;
        if (wideWordIndex >= wideTouchedAtomWords.size())
        {
            wideTouchedAtomWords.resize(wideWordIndex + 1, 0);
        }
        uint64_t &word = wideTouchedAtomWords[wideWordIndex];
        if (word == 0) touchedWideAtomWordIndices.push_back(wordIndex);
        word |= bit;
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
     * @param edgeCounts Edge counts parallel to the appended output masks
     * @param tempMaskList Reusable component-mask buffer; must not alias maskList
     */
    [[gnu::noinline]] void splitSmallWithBuffers(
        vector<EdgeMask> &maskList,
        vi &edgeCounts,
        vector<EdgeMask> &tempMaskList
    )
    {
        static_cast<void>(tempMaskList);
        componentMaskWords.clear();
        auto addTouchedAtom = [&](size_t index) {
            find(index);
            const int root = elements[index].parent;
            int32_t &component = elements[root].component;
            if (component == -1)
            {
                component = static_cast<int32_t>(componentMaskWords.size());
                componentMaskWords.push_back(0);
            }
            componentMaskWords[component] |=
                uint64_t{1} << elements[index].val;
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
        if (!touchedWideAtomWordIndices.empty()) [[unlikely]]
        {
            if (!is_sorted(
                touchedWideAtomWordIndices.begin(),
                touchedWideAtomWordIndices.end()
            ))
            {
                sort(
                    touchedWideAtomWordIndices.begin(),
                    touchedWideAtomWordIndices.end()
                );
            }
            for (const size_t wordIndex : touchedWideAtomWordIndices)
            {
                atoms = wideTouchedAtomWords[wordIndex - inlineAtomWordCount];
                while (atoms != 0)
                {
                    const size_t bitIndex = std::countr_zero(atoms);
                    addTouchedAtom(wordIndex * atomWordBits + bitIndex);
                    atoms &= atoms - 1;
                }
            }
        }
        for (const pii &extra : extraVals)
        {
            const size_t root = find(extra.first);
            componentMaskWords[elements[root].component] |=
                uint64_t{1} << extra.second;
        }
        for (const uint64_t word : componentMaskWords)
        {
            const int edgeCount = static_cast<int>(std::popcount(word));
            if (edgeCount > 1)
            {
                maskList.emplace_back(word);
                edgeCounts.push_back(edgeCount);
            }
        }
    }

    [[gnu::always_inline]] void splitWithBuffers(
        vector<EdgeMask> &maskList,
        vi &edgeCounts,
        vector<EdgeMask> &tempMaskList
    )
    {
        if (EdgeMask::activeWordCount() <= 1) [[likely]]
        {
            splitSmallWithBuffers(maskList, edgeCounts, tempMaskList);
            return;
        }
        if (EdgeMask::activeWordCount() == 2) [[likely]]
        {
            splitTwoWordWithBuffers(maskList, edgeCounts, tempMaskList);
            return;
        }
        splitWideWithBuffers(maskList, edgeCounts, tempMaskList);
    }

    [[gnu::noinline]] void splitTwoWordWithBuffers(
        vector<EdgeMask> &maskList,
        vi &edgeCounts,
        vector<EdgeMask> &tempMaskList
    )
    {
        static_cast<void>(tempMaskList);
        componentMaskWords.clear();
        auto setComponentEdge = [&](size_t componentOffset, size_t edge) {
            componentMaskWords[componentOffset + edge / 64] |=
                uint64_t{1} << (edge % 64);
        };
        auto addTouchedAtom = [&](size_t index) {
            find(index);
            const int root = elements[index].parent;
            int32_t &component = elements[root].component;
            if (component == -1)
            {
                component = static_cast<int32_t>(componentMaskWords.size());
                componentMaskWords.push_back(0);
                componentMaskWords.push_back(0);
            }
            setComponentEdge(component, elements[index].val);
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
        if (!touchedWideAtomWordIndices.empty()) [[unlikely]]
        {
            if (!is_sorted(
                touchedWideAtomWordIndices.begin(),
                touchedWideAtomWordIndices.end()
            ))
            {
                sort(
                    touchedWideAtomWordIndices.begin(),
                    touchedWideAtomWordIndices.end()
                );
            }
            for (const size_t wordIndex : touchedWideAtomWordIndices)
            {
                atoms = wideTouchedAtomWords[wordIndex - inlineAtomWordCount];
                while (atoms != 0)
                {
                    const size_t bitIndex = std::countr_zero(atoms);
                    addTouchedAtom(wordIndex * atomWordBits + bitIndex);
                    atoms &= atoms - 1;
                }
            }
        }
        for (const pii &extra : extraVals)
        {
            const size_t root = find(extra.first);
            setComponentEdge(elements[root].component, extra.second);
        }
        for (size_t offset = 0;
             offset < componentMaskWords.size();
             offset += 2)
        {
            const uint64_t *words = componentMaskWords.data() + offset;
            const int edgeCount = static_cast<int>(
                std::popcount(words[0]) + std::popcount(words[1])
            );
            if (edgeCount > 1)
            {
                maskList.push_back(EdgeMask::fromActiveWords(words));
                edgeCounts.push_back(edgeCount);
            }
        }
    }

    [[gnu::noinline]] void splitWideWithBuffers(
        vector<EdgeMask> &maskList,
        vi &edgeCounts,
        vector<EdgeMask> &tempMaskList
    )
    {
        tempMaskList.clear();
        componentMaskWords.clear();
        const size_t edgeWordCount = EdgeMask::activeWordCount();
        auto setComponentEdge = [&](size_t componentOffset, size_t edge) {
            componentMaskWords[componentOffset + edge / 64] |=
                uint64_t{1} << (edge % 64);
        };
        auto addTouchedAtom = [&](size_t index) {
            find(index);
            const int root = elements[index].parent;
            int32_t &component = elements[root].component;
            if (component == -1)
            {
                if (componentMaskWords.size() >
                    static_cast<size_t>(numeric_limits<int32_t>::max()))
                {
                    throw length_error("too many residual component words");
                }
                component = static_cast<int32_t>(componentMaskWords.size());
                componentMaskWords.resize(
                    componentMaskWords.size() + edgeWordCount,
                    0
                );
            }
            setComponentEdge(component, elements[index].val);
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
        if (!touchedWideAtomWordIndices.empty()) [[unlikely]]
        {
            if (!is_sorted(
                touchedWideAtomWordIndices.begin(),
                touchedWideAtomWordIndices.end()
            ))
            {
                sort(
                    touchedWideAtomWordIndices.begin(),
                    touchedWideAtomWordIndices.end()
                );
            }
            for (const size_t wordIndex : touchedWideAtomWordIndices)
            {
                atoms = wideTouchedAtomWords[wordIndex - inlineAtomWordCount];
                while (atoms != 0)
                {
                    const size_t bitIndex = std::countr_zero(atoms);
                    addTouchedAtom(wordIndex * atomWordBits + bitIndex);
                    atoms &= atoms - 1;
                }
            }
        }
        for (size_t i = 0; i < extraVals.size(); i++)
        {
            const size_t root = find(extraVals[i].first);
            const int32_t component = elements[root].component;
            setComponentEdge(component, extraVals[i].second);
        }
        for (size_t offset = 0;
             offset < componentMaskWords.size();
             offset += edgeWordCount)
        {
            const uint64_t *words = componentMaskWords.data() + offset;
            int edgeCount = 0;
            for (size_t wordIndex = 0; wordIndex < edgeWordCount; wordIndex++)
            {
                edgeCount += static_cast<int>(std::popcount(words[wordIndex]));
            }
            if (edgeCount > 1)
            {
                maskList.push_back(EdgeMask::fromActiveWords(words));
                edgeCounts.push_back(edgeCount);
            }
        }
    }
};
