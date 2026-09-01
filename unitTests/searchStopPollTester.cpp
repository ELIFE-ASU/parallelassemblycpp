#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <csignal>
#include <ctime>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
using IntegerVector = vector<int>;
using BooleanVector = vector<bool>;
using IntegerPair = pair<int, int>;

#include "../v5/activeWordMask.h"
#include "../v5/distributedRootMapping.h"
#include "../v5/globalPrimitives.h"

bool distributedRootMappingHasExactCoverage(
    size_t rootJobCount,
    size_t leaseSize,
    size_t stripeCount
)
{
    if (leaseSize == 0 || stripeCount == 0) return false;
    if (stripeCount > numeric_limits<size_t>::max() / leaseSize) return false;
    const size_t blockWidth = leaseSize * stripeCount;
    const size_t remainder = rootJobCount % blockWidth;
    const size_t padding = remainder == 0 ? 0 : blockWidth - remainder;
    if (rootJobCount > numeric_limits<size_t>::max() - padding) return false;
    const size_t paddedSlotCount = rootJobCount + padding;

    vector<unsigned char> visits(rootJobCount, 0);
    size_t paddingSlots = 0;
    for (size_t ordinal = 0; ordinal < paddedSlotCount; ++ordinal)
    {
        const size_t rootJobIndex = distributedStripedRootJobIndex(
            ordinal,
            leaseSize,
            stripeCount,
            rootJobCount
        );
        if (rootJobIndex >= rootJobCount)
        {
            ++paddingSlots;
            continue;
        }
        if (visits[rootJobIndex] != 0) return false;
        visits[rootJobIndex] = 1;
    }
    return
        paddingSlots == padding &&
        accumulate(visits.begin(), visits.end(), size_t{0}) == rootJobCount;
}

bool distributedRootMappingsAreValid()
{
    // Exhaust the compact layouts most likely to expose partial final chunks.
    for (size_t rootJobCount = 0; rootJobCount <= 65; ++rootJobCount)
    {
        for (size_t leaseSize = 1; leaseSize <= 9; ++leaseSize)
        {
            for (size_t stripeCount = 1; stripeCount <= 8; ++stripeCount)
            {
                if (!distributedRootMappingHasExactCoverage(
                    rootJobCount,
                    leaseSize,
                    stripeCount
                )) return false;
            }
        }
    }

    // Production-shaped odd and padded tails, including fixed leases > 1.
    constexpr array<array<size_t, 3>, 8> cases{{
        {{9, 7, 4}},
        {{17, 3, 5}},
        {{45, 2, 3}},
        {{45, 7, 2}},
        {{45, 7, 3}},
        {{45, 7, 4}},
        {{97, 8, 7}},
        {{238, 7, 4}}
    }};
    for (const auto &[rootJobCount, leaseSize, stripeCount] : cases)
    {
        if (!distributedRootMappingHasExactCoverage(
            rootJobCount,
            leaseSize,
            stripeCount
        )) return false;
    }
    return true;
}

void setInterruptFlag(bool value)
{
    searchCancellationFlag.store(value);
#ifdef _WIN32
    interruptFlag.store(value);
#else
    interruptFlag = value;
#endif
}

int main()
{
    if (!distributedRootMappingsAreValid()) return 1;

    maximumRuntimeTicks = numeric_limits<unsigned long long>::max();
    searchStopPollCountdown = 0;
    searchStopInnerPollCountdown = 0;
    setInterruptFlag(false);

    // The first inner-loop call polls immediately.
    if (searchShouldStopPeriodically()) return 1;

    // A newly requested stop may be skipped only until the next 128th poll.
    setInterruptFlag(true);
    for (size_t skipped = 1; skipped < searchStopPollInterval; skipped++)
    {
        if (searchShouldStopPeriodically()) return 1;
    }
    if (!searchShouldStopPeriodically()) return 1;

    // Once observed, a stop remains immediate while callers unwind.
    if (!searchShouldStopPeriodically()) return 1;
    if (searchStopInnerPollCountdown != 0) return 1;

    // Search boundaries continue to observe interrupts immediately.
    searchStopPollCountdown = searchStopPollInterval - 1;
    if (!searchShouldStop()) return 1;

    // Ordinary search boundaries arm and poll at the full production cadence.
    setInterruptFlag(false);
    runtimeLimitReached = false;
    startTime = clock();
    maximumRuntimeTicks = numeric_limits<unsigned long long>::max() - 1;
    searchStopPollCountdown = 0;
    if (searchShouldStop()) return 1;
    if (searchStopPollCountdown != searchStopPollInterval - 1) return 1;

    maximumRuntimeTicks = 0;
    for (size_t skipped = 1; skipped < searchStopPollInterval; skipped++)
    {
        if (searchShouldStop()) return 1;
    }
    if (!searchShouldStop()) return 1;
    if (!runtimeLimitReached) return 1;

    // A due inner-loop poll forces a deadline sample.
    setInterruptFlag(false);
    runtimeLimitReached = false;
    searchStopPollCountdown = searchStopPollInterval - 1;
    searchStopInnerPollCountdown = 2;
    if (searchShouldStopPeriodically()) return 1;
    if (searchShouldStopPeriodically()) return 1;
    if (!searchShouldStopPeriodically()) return 1;
    if (!runtimeLimitReached) return 1;

    // A direct zero-runtime check retains its immediate-stop semantics.
    setInterruptFlag(false);
    runtimeLimitReached = false;
    searchStopPollCountdown = 0;
    if (!searchShouldStop()) return 1;
    if (!runtimeLimitReached) return 1;

    return 0;
}
