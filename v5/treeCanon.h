using treeCanonNodeId = std::uint64_t;
using treeCanonAtomId = std::uint64_t;

/**
 * @brief Canonical identity of an unrooted, labelled tree.
 *
 * A tree has either one centroid, represented by first alone, or two
 * centroids, represented by the sorted first/second pair and their connecting
 * bond. Zero is reserved for graphs without a tree canonical form.
 */
struct treeCanonForm
{
    treeCanonNodeId first = 0;
    treeCanonNodeId second = 0;
    unsigned char centralBond = 0;

    bool empty() const {return first == 0;}
    bool operator==(const treeCanonForm &) const = default;
};

/**
 * @brief One labelled edge to an already-canonical child subtree.
 */
struct treeCanonChild
{
    unsigned char bondType;
    treeCanonNodeId subtree;

    bool operator==(const treeCanonChild &) const = default;
};

/**
 * @brief Exact structural key for a rooted subtree.
 */
struct treeCanonSignature
{
    treeCanonAtomId atomType;
    std::vector<treeCanonChild> children;

    bool operator==(const treeCanonSignature &) const = default;
};

struct treeCanonSignatureHash
{
    static void combine(std::size_t &seed, std::size_t value)
    {
        seed ^= value + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL)
            + (seed << 6) + (seed >> 2);
    }

    std::size_t operator()(const treeCanonSignature &signature) const
    {
        std::size_t result = std::hash<treeCanonAtomId>{}(signature.atomType);
        for (const treeCanonChild &child : signature.children)
        {
            combine(result, std::hash<unsigned char>{}(child.bondType));
            combine(result, std::hash<treeCanonNodeId>{}(child.subtree));
        }
        return result;
    }
};

/**
 * IDs are opaque and stable within one interner generation. Production
 * discards graphHashMap before resetting these interners between calculations;
 * forms from different generations must not be mixed.
 */
inline std::unordered_map<std::string, treeCanonAtomId> treeCanonAtomInterner;
inline std::vector<treeCanonNodeId> treeCanonLeafInterner;
inline std::unordered_map<
    treeCanonSignature,
    treeCanonNodeId,
    treeCanonSignatureHash
> treeCanonInterner;

void clearTreeCanonInterner()
{
    treeCanonInterner.clear();
    treeCanonAtomInterner.clear();
    treeCanonLeafInterner.clear();
}

treeCanonAtomId internTreeCanonAtom(const std::string &atomType)
{
    const treeCanonAtomId nextId =
        static_cast<treeCanonAtomId>(treeCanonAtomInterner.size()) + 1;
    return treeCanonAtomInterner.try_emplace(atomType, nextId).first->second;
}

treeCanonNodeId internTreeCanonSignature(treeCanonSignature signature)
{
    const treeCanonNodeId nextId =
        static_cast<treeCanonNodeId>(treeCanonInterner.size()) + 1;
    return treeCanonInterner.try_emplace(
        std::move(signature),
        nextId
    ).first->second;
}

treeCanonNodeId internTreeCanonNode(
    treeCanonAtomId atomType,
    std::vector<treeCanonChild> children
)
{
    std::sort(
        children.begin(),
        children.end(),
        [](const treeCanonChild &left, const treeCanonChild &right)
        {
            if (left.bondType != right.bondType)
                return left.bondType < right.bondType;
            return left.subtree < right.subtree;
        }
    );
    if (children.empty())
    {
        if (treeCanonLeafInterner.size() <= atomType)
            treeCanonLeafInterner.resize(atomType + 1, 0);
        treeCanonNodeId &leaf = treeCanonLeafInterner[atomType];
        if (leaf == 0)
            leaf = internTreeCanonSignature({atomType, {}});
        return leaf;
    }
    return internTreeCanonSignature({atomType, std::move(children)});
}

treeCanonNodeId internTreeCanonNode(
    molGraph &mg,
    int node,
    std::vector<treeCanonChild> children
)
{
    return internTreeCanonNode(
        internTreeCanonAtom(mg.mg[node].type),
        std::move(children)
    );
}

/**
 * @brief Iteratively canonicalise an acyclic molecular graph.
 *
 * AHU leaf peeling processes a complete layer at a time. Each removed node is
 * represented by an exact interned ID built from its atom label and sorted
 * (bond-byte, child-ID) multiset. The one or two unremoved nodes are the tree
 * centroids, so no recursive traversal or subtree-string concatenation is
 * required.
 *
 * @param mg Acyclic, connected molGraph
 * @param n Retained for API compatibility; centroid selection is root-free.
 */
treeCanonForm centroidTreeCanon(molGraph &mg, int n)
{
    (void)n;
    const std::size_t nodeCount = mg.mg.size();
    if (nodeCount == 0) return {};
    if (
        mg.totalBonds < 0 ||
        static_cast<std::size_t>(mg.totalBonds) != nodeCount - 1
    ) return {};

    std::vector<int> remainingDegree(nodeCount, 0);
    std::vector<unsigned char> removed(nodeCount, 0);
    std::vector<std::vector<treeCanonChild>> children(nodeCount);
    std::vector<int> leaves;
    leaves.reserve(nodeCount);
    std::size_t remaining = 0;

    for (std::size_t node = 0; node < nodeCount; node++)
    {
        // X was a legacy deletion sentinel with inconsistent partial-tree
        // semantics. Route such input through exact whole-graph canonicalisation.
        if (mg.mg[node].type == "X") return {};
        remaining++;
        remainingDegree[node] = static_cast<int>(mg.mg[node].list.size());
        if (remainingDegree[node] <= 1)
            leaves.push_back(static_cast<int>(node));
    }

    while (remaining > 2)
    {
        // A connected acyclic graph always has a leaf. Returning no tree form
        // safely routes malformed input through exact whole-graph canonicalisation.
        if (leaves.empty()) return {};
        remaining -= leaves.size();
        for (const int leaf : leaves) removed[leaf] = 1;

        std::vector<int> nextLeaves;
        for (const int leaf : leaves)
        {
            const treeCanonNodeId subtree = internTreeCanonNode(
                mg,
                leaf,
                std::move(children[leaf])
            );
            for (const bond &edge : mg.mg[leaf].list)
            {
                const int neighbour = edge.n;
                if (
                    neighbour < 0 ||
                    static_cast<std::size_t>(neighbour) >= nodeCount ||
                    removed[neighbour]
                ) continue;

                children[neighbour].push_back({
                    static_cast<unsigned char>(static_cast<char>(edge.type)),
                    subtree
                });
                remainingDegree[neighbour]--;
                if (remainingDegree[neighbour] == 1)
                    nextLeaves.push_back(neighbour);
                break;
            }
        }
        leaves = std::move(nextLeaves);
    }

    std::vector<int> centroids;
    centroids.reserve(2);
    for (std::size_t node = 0; node < nodeCount; node++)
    {
        if (!removed[node]) centroids.push_back(static_cast<int>(node));
    }
    if (centroids.empty() || centroids.size() > 2) return {};

    const treeCanonNodeId first = internTreeCanonNode(
        mg,
        centroids[0],
        std::move(children[centroids[0]])
    );
    if (centroids.size() == 1) return {first, 0, 0};

    unsigned char centralBond = 0;
    bool centralBondFound = false;
    for (const bond &edge : mg.mg[centroids[0]].list)
    {
        if (edge.n == centroids[1])
        {
            centralBond = static_cast<unsigned char>(static_cast<char>(edge.type));
            centralBondFound = true;
            break;
        }
    }
    if (!centralBondFound) return {};

    const treeCanonNodeId second = internTreeCanonNode(
        mg,
        centroids[1],
        std::move(children[centroids[1]])
    );
    if (first <= second) return {first, second, centralBond};
    return {second, first, centralBond};
}
