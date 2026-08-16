#include <cstdint>
#include <limits>
#include <type_traits>
#ifdef _WIN32
    #include <atomic>
#endif

template <typename T1, typename T2, typename T3>
struct triple
{
    T1 a; T2 b; T3 c;
    triple(){}
    triple(T1 &_a, T2 &_b, T3 &_c): a(_a), b(_b), c(_c) {}
};

int DISCHARGE_FREQUENCY = 30000000, ENUM_MAX = 50000000;
int minAIfound = -1;
int recursiveCount = 1;
int assemblyIx = 0;

std::unordered_map<string, int> atypeHash;
vector<double> coords;
string moleculeName;
standardBitset allEdges;
#ifdef _WIN32
    std::atomic_bool interruptFlag = false;
    std::atomic_bool userInterruptReceived = false;
#else
    volatile std::sig_atomic_t interruptFlag = 0;
    volatile std::sig_atomic_t userInterruptReceived = 0;
#endif
clock_t startTime = 0;
unsigned long long runTimeMax = 18446744073709551615ULL;
bool runtimeLimitReached = false;
bool enumerationLimitReached = false;

typedef triple<int, int, int> iii;
typedef triple<short, short, short> edgeL;
unsigned int totalBonds = 0;
vector<edgeL> removedEdges;
vector<edgeL> originalEdgeList, univEdgeList;

/// Hash table for edgelists for pathway algorithm
std::unordered_map<standardBitset, pii> bitsetHashTable;

bool isPathway = 1, removeHydrogens = 1, disjointCompensation = 0, memTest = 0, writeIntermediateMAs = 0;
int disjointFragments = 1;
vector<pair<unsigned long long, int> > intermediateMAs;

bool interruptionRequested()
{
    #ifdef _WIN32
        return interruptFlag.load();
    #else
        return interruptFlag != 0;
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
bool searchShouldStop()
{
    if (interruptionRequested()) return true;
    if (runTimeMax == std::numeric_limits<unsigned long long>::max()) return false;
    if (elapsedClockTicks() < runTimeMax) return false;

    runtimeLimitReached = true;
    #ifdef _WIN32
        interruptFlag.store(true);
    #else
        interruptFlag = 1;
    #endif
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
