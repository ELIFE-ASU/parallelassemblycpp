#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <csignal>
#include <ctime>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
using vi = vector<int>;
using vb = vector<bool>;
using pii = pair<int, int>;

#include "../v5/activeWordMask.h"
#include "../v5/globalPrimitives.h"

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
    runTimeMax = numeric_limits<unsigned long long>::max();
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
    runTimeMax = numeric_limits<unsigned long long>::max() - 1;
    searchStopPollCountdown = 0;
    if (searchShouldStop()) return 1;
    if (searchStopPollCountdown != searchStopPollInterval - 1) return 1;

    runTimeMax = 0;
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
