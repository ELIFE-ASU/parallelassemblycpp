#pragma once

#include <ctime>
#include <limits>
#include <type_traits>

namespace assembly_clock
{
using UnsignedClock = std::make_unsigned_t<std::clock_t>;

/** Subtract clock samples modulo clock_t's unsigned range; -1 means unavailable. */
[[nodiscard]] constexpr UnsignedClock difference(
    std::clock_t start,
    std::clock_t end
) noexcept
{
    constexpr auto unavailable = static_cast<std::clock_t>(-1);
    if (start == unavailable || end == unavailable) return 0;
    return static_cast<UnsignedClock>(end) - static_cast<UnsignedClock>(start);
}

/** Runtime budgets use unsigned long long even on platforms with wider clocks. */
[[nodiscard]] constexpr unsigned long long budgetTicks(
    std::clock_t start,
    std::clock_t end
) noexcept
{
    const UnsignedClock elapsed = difference(start, end);
    if constexpr (sizeof(UnsignedClock) > sizeof(unsigned long long))
    {
        if (elapsed > static_cast<UnsignedClock>(
            std::numeric_limits<unsigned long long>::max()
        )) return std::numeric_limits<unsigned long long>::max();
    }
    return static_cast<unsigned long long>(elapsed);
}
}
