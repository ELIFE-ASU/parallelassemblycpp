#include <atomic>

template <typename T1, typename T2, typename T3>
struct triple
{
    T1 a; T2 b; T3 c;
    triple(T1 &_a, T2 &_b, T3 &_c): a(_a), b(_b), c(_c) {}
};

constexpr int unknownCanonicalId = -1;

/** @brief One assembly fragment and the metadata proved for its exact mask. */
struct assemblyFragment
{
    EdgeMask mask;
    int32_t canonicalId = unknownCanonicalId;
    uint32_t edgeCount : 31 = 0;
    uint32_t connected : 1 = false;

    assemblyFragment() = default;

    assemblyFragment(
        const EdgeMask &_mask,
        int _edgeCount,
        int _canonicalId = unknownCanonicalId,
        bool _connected = false
    ):
        mask(_mask),
        canonicalId(_canonicalId),
        edgeCount(static_cast<uint32_t>(_edgeCount)),
        connected(_connected) {}

    assemblyFragment(
        EdgeMask &&_mask,
        int _edgeCount,
        int _canonicalId = unknownCanonicalId,
        bool _connected = false
    ):
        mask(std::move(_mask)),
        canonicalId(_canonicalId),
        edgeCount(static_cast<uint32_t>(_edgeCount)),
        connected(_connected) {}

    /** Retain allowed edges, invalidating metadata only when the mask changes. */
    [[gnu::always_inline]] bool retainEdges(const EdgeMask &allowedMask)
    {
        if (allowedMask.contains(mask)) return false;
        mask &= allowedMask;
        edgeCount = static_cast<uint32_t>(mask.count());
        canonicalId = unknownCanonicalId;
        connected = false;
        return true;
    }
};

static_assert(sizeof(assemblyFragment) == 2 * sizeof(uint64_t));

int ENUM_MAX = 50000000;

ASSEMBLYCPP_SEARCH_LOCAL string moleculeName;
ASSEMBLYCPP_SEARCH_LOCAL EdgeMask allEdges;
#ifdef _WIN32
    std::atomic_bool interruptFlag = false;
    std::atomic_bool userInterruptReceived = false;
#else
    volatile std::sig_atomic_t interruptFlag = 0;
    volatile std::sig_atomic_t userInterruptReceived = 0;
#endif
// Solver workers use a real C++ atomic for cooperative cancellation. POSIX
// signal handlers retain sig_atomic_t flags and never write this object.
std::atomic_bool searchCancellationFlag = false;
ASSEMBLYCPP_SEARCH_LOCAL clock_t startTime = 0;
unsigned long long runTimeMax = std::numeric_limits<unsigned long long>::max();
ASSEMBLYCPP_SEARCH_LOCAL bool runtimeLimitReached = false;
ASSEMBLYCPP_SEARCH_LOCAL bool enumerationLimitReached = false;

using edgeL = triple<short, short, short>;
ASSEMBLYCPP_SEARCH_LOCAL unsigned int totalBonds = 0;
ASSEMBLYCPP_SEARCH_LOCAL vector<edgeL> originalEdgeList, univEdgeList;

/// Hash table for edgelists for pathway algorithm
ASSEMBLYCPP_SEARCH_LOCAL std::unordered_map<EdgeMask, pii> bitsetHashTable;

bool isPathway = true;
bool removeHydrogens = true;
bool verbose = false;
bool disjointCompensation = false;
bool memTest = false;
bool writeIntermediateMAs = false;
ASSEMBLYCPP_SEARCH_LOCAL int lastCalculatedAssemblyIndex = -1;
ASSEMBLYCPP_SEARCH_LOCAL int disjointFragments = 1;
ASSEMBLYCPP_SEARCH_LOCAL vector<pair<unsigned long long, int>> intermediateMAs;

// Parallel builds replicate the complete solver state. MPI ranks retain a
// deterministic modulo partition of root branches, while workers inside a
// rank dynamically claim rank-local ranges from one shared cursor.
ASSEMBLYCPP_SEARCH_LOCAL size_t searchRankPartitionIndex = 0;
ASSEMBLYCPP_SEARCH_LOCAL size_t searchRankPartitionCount = 1;
ASSEMBLYCPP_SEARCH_LOCAL size_t searchRootBranchOrdinal = 0;
ASSEMBLYCPP_SEARCH_LOCAL size_t searchBranchLeaseSize = 1;
ASSEMBLYCPP_SEARCH_LOCAL std::atomic<size_t> *sharedBranchLeaseCursor = nullptr;
ASSEMBLYCPP_SEARCH_LOCAL size_t searchBranchLeaseCount = 0;
ASSEMBLYCPP_SEARCH_LOCAL size_t searchBranchAssignmentCount = 0;
ASSEMBLYCPP_SEARCH_LOCAL std::atomic<int> *sharedAssemblyIndex = nullptr;
ASSEMBLYCPP_SEARCH_LOCAL bool suppressSearchOutput = false;

bool interruptionRequested()
{
    #ifdef _WIN32
        return interruptFlag.load() || searchCancellationFlag.load();
    #else
        return interruptFlag != 0 || searchCancellationFlag.load();
    #endif
}

bool receivedUserInterrupt()
{
    #ifdef _WIN32
        return userInterruptReceived.load();
    #else
        return userInterruptReceived != 0;
    #endif
}

/**
 * @brief Return elapsed std::clock ticks without signed arithmetic.
 */
unsigned long long elapsedClockTicks()
{
    const clock_t now = clock();
    const clock_t clockError = static_cast<clock_t>(-1);
    if (now == clockError || startTime == clockError) return 0;

    using unsignedClock = std::make_unsigned_t<clock_t>;
    const unsignedClock elapsed =
        static_cast<unsignedClock>(now) - static_cast<unsignedClock>(startTime);
    const std::uintmax_t elapsedWide = static_cast<std::uintmax_t>(elapsed);
    const std::uintmax_t outputMax =
        static_cast<std::uintmax_t>(std::numeric_limits<unsigned long long>::max());
    if (elapsedWide > outputMax)
        return std::numeric_limits<unsigned long long>::max();
    return static_cast<unsigned long long>(elapsedWide);
}

/**
 * @brief Cooperatively stop the search once its std::clock budget is spent.
 */
constexpr size_t searchStopPollInterval = 128;
static_assert(std::has_single_bit(searchStopPollInterval));
ASSEMBLYCPP_SEARCH_LOCAL size_t searchStopPollCountdown = 0;
ASSEMBLYCPP_SEARCH_LOCAL size_t searchStopInnerPollCountdown = 0;

bool searchShouldStop()
{
    if (interruptionRequested()) return true;
    if (runTimeMax == std::numeric_limits<unsigned long long>::max()) return false;

    if (searchStopPollCountdown != 0)
    {
        --searchStopPollCountdown;
        return false;
    }

    searchStopPollCountdown = searchStopPollInterval - 1;
    if (elapsedClockTicks() < runTimeMax) return false;

    runtimeLimitReached = true;
    searchStopPollCountdown = 0;
    searchCancellationFlag.store(true);
    return true;
}

/**
 * @brief Check cancellation at a bounded cadence inside cheap inner loops.
 */
bool searchShouldStopPeriodically()
{
    if (searchStopInnerPollCountdown != 0)
    {
        --searchStopInnerPollCountdown;
        return false;
    }

    searchStopInnerPollCountdown = searchStopPollInterval - 1;
    // A due inner-loop poll must sample the deadline even if ordinary search
    // boundaries recently did so.
    searchStopPollCountdown = 0;
    if (!searchShouldStop()) return false;

    searchStopInnerPollCountdown = 0;
    return true;
}

/**
 * @brief Number of extra disconnected components used for index compensation.
 */
int disjointComponentAdjustment()
{
    return disjointFragments > 1 ? disjointFragments - 1 : 0;
}

int compensateDisjointAssemblyIndex(int index)
{
    return disjointCompensation ? index - disjointComponentAdjustment() : index;
}
