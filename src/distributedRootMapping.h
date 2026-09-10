#pragma once

#include <cstddef>
#include <limits>

/**
 * Map one padded distributed-queue slot to its striped root-job index.
 *
 * Queue chunks are contiguous so they can be leased cheaply. Transposing the
 * chunks across stripes preserves the former small-stride search order while
 * allowing any rank to request the next chunk. Indices at or above the real
 * root count are padding and must be skipped by the caller.
 */
[[nodiscard]] constexpr std::size_t distributedStripedRootJobIndex(
    std::size_t queueOrdinal,
    std::size_t leaseSize,
    std::size_t stripeCount,
    std::size_t rootJobCount
) noexcept
{
    if (leaseSize == 0 || stripeCount == 0) return rootJobCount;

    const std::size_t chunk = queueOrdinal / leaseSize;
    const std::size_t offset = queueOrdinal % leaseSize;
    const std::size_t stripe = chunk % stripeCount;
    const std::size_t block = chunk / stripeCount;
    if (
        block >
        (std::numeric_limits<std::size_t>::max() - offset) / leaseSize
    ) return rootJobCount;
    const std::size_t stripeOrdinal = block * leaseSize + offset;
    if (
        stripeOrdinal >
        (std::numeric_limits<std::size_t>::max() - stripe) / stripeCount
    ) return rootJobCount;
    return stripeOrdinal * stripeCount + stripe;
}
