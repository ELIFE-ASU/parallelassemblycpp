#include <cassert>
#include <ctime>
#include <limits>

#include "../src/clockTicks.h"

int main()
{
    constexpr auto unavailable = static_cast<std::clock_t>(-1);
    static_assert(assembly_clock::difference(10, 25) == 15);
    static_assert(assembly_clock::difference(25, 25) == 0);
    static_assert(assembly_clock::difference(unavailable, 25) == 0);
    static_assert(assembly_clock::difference(25, unavailable) == 0);
    static_assert(assembly_clock::budgetTicks(10, 25) == 15);
    static_assert(assembly_clock::budgetTicks(unavailable, unavailable) == 0);

    // Exercise rollover without signed subtraction or relying on a running clock.
    constexpr auto before = std::numeric_limits<std::clock_t>::max() - 2;
    constexpr auto after = std::numeric_limits<std::clock_t>::lowest() + 2;
    static_assert(assembly_clock::difference(before, after) == 5);
    assert(assembly_clock::budgetTicks(before, after) == 5);
}
