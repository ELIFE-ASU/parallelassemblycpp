#ifndef ASSEMBLY_TRANSPOSITION_TABLE_H
#define ASSEMBLY_TRANSPOSITION_TABLE_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <memory_resource>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

/**
 * @brief Compact score-only transposition table for canonical assembly keys.
 *
 * Lookups accept borrowed spans. A key is copied into monotonic arena storage
 * only when it is first inserted; hits neither allocate nor retain the span.
 * The open-addressed index has no erase operation, so it needs no tombstones.
 */
class assemblyTranspositionTable
{
public:
    enum class result
    {
        inserted,
        improved,
        dominated
    };

    explicit assemblyTranspositionTable(
        std::size_t initialCapacity = minimumCapacity,
        std::pmr::memory_resource *keyUpstream =
            std::pmr::get_default_resource()
    ):
        keyArena(requireUpstream(keyUpstream)),
        slots(normaliseCapacity(initialCapacity))
    {}

    assemblyTranspositionTable(const assemblyTranspositionTable &) = delete;
    assemblyTranspositionTable &operator=(
        const assemblyTranspositionTable &
    ) = delete;
    assemblyTranspositionTable(assemblyTranspositionTable &&) = delete;
    assemblyTranspositionTable &operator=(
        assemblyTranspositionTable &&
    ) = delete;

    /**
     * @brief Insert or improve the duplicated-bond score for a canonical key.
     *
     * @return inserted for a new key, improved for a strictly larger score,
     * or dominated when an equal or larger score was already stored.
     */
    result consider(std::span<const int> key, int sumDupBonds)
    {
        validateKeyLength(key.size());
        const std::uint32_t hash = hashKey(key);
        const findResult found = find(key, hash);
        if (found.distance == foundDistance)
        {
            slot *existing = &slots[found.index];
            if (sumDupBonds <= existing->bestSumDupBonds)
                return result::dominated;
            existing->bestSumDupBonds = sumDupBonds;
            return result::improved;
        }

        return insertMiss(key, hash, sumDupBonds, found);
    }

    [[nodiscard]] std::size_t size() const noexcept
    {
        return sizeValue;
    }

    [[nodiscard]] std::size_t capacity() const noexcept
    {
        return slots.size();
    }

private:
    static constexpr std::size_t minimumCapacity = 8;

    struct alignas(int) storedKeyHeader
    {
        std::uint32_t length;
    };

    struct alignas(16) slot
    {
        const storedKeyHeader *key = nullptr;
        std::uint32_t hash = 0;
        int bestSumDupBonds = 0;
    };

    static_assert(sizeof(slot) == 16);
    static_assert(std::is_trivially_destructible_v<slot>);

    struct findResult
    {
        std::size_t index;
        std::size_t distance;
    };

    static constexpr std::size_t foundDistance =
        std::numeric_limits<std::size_t>::max();

    std::pmr::monotonic_buffer_resource keyArena;
    std::vector<slot> slots;
    std::size_t sizeValue = 0;

    static std::pmr::memory_resource *requireUpstream(
        std::pmr::memory_resource *resource
    )
    {
        if (resource == nullptr)
            throw std::invalid_argument("key arena upstream is null");
        return resource;
    }

    static std::size_t normaliseCapacity(std::size_t requested)
    {
        std::size_t result = minimumCapacity;
        const std::size_t maximumCapacity =
            std::numeric_limits<std::size_t>::max() / sizeof(slot);
        while (result < requested)
        {
            if (result > maximumCapacity / 2)
                throw std::length_error("transposition table is too large");
            result *= 2;
        }
        return result;
    }

    static void validateKeyLength(std::size_t length)
    {
        if (length > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("assembly-state key is too long");
    }

    static std::size_t maximumEntries(std::size_t capacity) noexcept
    {
        // Keep at least one fifth of the slots empty. The subtraction form
        // avoids overflow for large capacities.
        return capacity - (capacity + 4) / 5;
    }

    static std::uint32_t hashKey(std::span<const int> key) noexcept
    {
        std::uint32_t hash =
            0x811c9dc5U ^ static_cast<std::uint32_t>(key.size());
        for (const int value : key)
        {
            hash ^= static_cast<std::uint32_t>(value);
            hash *= 0x01000193U;
        }
        // Avalanche the low bits used by the power-of-two index. This mix is
        // cheaper than a 64-bit finalizer for the short assembly-state keys.
        hash ^= hash >> 16;
        hash *= 0x7feb352dU;
        hash ^= hash >> 15;
        hash ^= hash >> 16;
        return hash;
    }

    static const int *keyValues(const storedKeyHeader *key) noexcept
    {
        return reinterpret_cast<const int *>(key + 1);
    }

    static bool keysEqual(
        const storedKeyHeader *stored,
        std::span<const int> candidate
    ) noexcept
    {
        return stored->length == candidate.size() && std::equal(
            candidate.begin(),
            candidate.end(),
            keyValues(stored)
        );
    }

    static std::size_t probeDistance(
        std::uint32_t hash,
        std::size_t index,
        std::size_t mask
    ) noexcept
    {
        return (index - (static_cast<std::size_t>(hash) & mask)) & mask;
    }

    findResult find(std::span<const int> key, std::uint32_t hash) noexcept
    {
        const std::size_t mask = slots.size() - 1;
        std::size_t index = static_cast<std::size_t>(hash) & mask;
        for (std::size_t distance = 0;; distance++)
        {
            slot &candidate = slots[index];
            if (candidate.key == nullptr)
                return {index, distance};
            if (
                candidate.hash == hash && keysEqual(candidate.key, key)
            ) return {index, foundDistance};
            // No resident entry can have a negative displacement, so avoid
            // calculating it on the overwhelmingly common first probe.
            if (
                distance != 0 &&
                probeDistance(candidate.hash, index, mask) < distance
            ) return {index, distance};
            index = (index + 1) & mask;
        }
    }

    const storedKeyHeader *copyKey(std::span<const int> key)
    {
        if (key.size_bytes() >
            std::numeric_limits<std::size_t>::max() - sizeof(storedKeyHeader))
        {
            throw std::length_error("assembly-state key allocation overflow");
        }
        const std::size_t bytes = sizeof(storedKeyHeader) + key.size_bytes();
        void *memory = keyArena.allocate(bytes, alignof(storedKeyHeader));
        storedKeyHeader *header = std::construct_at(
            static_cast<storedKeyHeader *>(memory),
            storedKeyHeader{static_cast<std::uint32_t>(key.size())}
        );
        std::uninitialized_copy(
            key.begin(),
            key.end(),
            reinterpret_cast<int *>(header + 1)
        );
        return header;
    }

    [[gnu::noinline]] result insertMiss(
        std::span<const int> key,
        std::uint32_t hash,
        int sumDupBonds,
        const findResult &found
    )
    {
        const bool needsGrowth =
            sizeValue >= maximumEntries(slots.size());
        if (needsGrowth) grow();

        const storedKeyHeader *storedKey = copyKey(key);
        slot incoming{storedKey, hash, sumDupBonds};
        if (needsGrowth) insertWithoutLookup(slots, incoming);
        else insertAt(slots, incoming, found.index, found.distance);
        ++sizeValue;
        return result::inserted;
    }

    static void insertWithoutLookup(
        std::vector<slot> &destination,
        slot incoming
    ) noexcept
    {
        const std::size_t mask = destination.size() - 1;
        std::size_t index = static_cast<std::size_t>(incoming.hash) & mask;
        std::size_t distance = 0;
        insertAt(destination, incoming, index, distance);
    }

    static void insertAt(
        std::vector<slot> &destination,
        slot incoming,
        std::size_t index,
        std::size_t distance
    ) noexcept
    {
        const std::size_t mask = destination.size() - 1;
        while (true)
        {
            slot &current = destination[index];
            if (current.key == nullptr)
            {
                current = incoming;
                return;
            }

            const std::size_t currentDistance =
                probeDistance(current.hash, index, mask);
            if (currentDistance < distance)
            {
                std::swap(current, incoming);
                distance = currentDistance;
            }
            index = (index + 1) & mask;
            ++distance;
        }
    }

    void grow()
    {
        if (slots.size() > slots.max_size() / 2)
            throw std::length_error("transposition table is too large");
        std::vector<slot> expanded(slots.size() * 2);
        for (const slot &entry : slots)
        {
            if (entry.key != nullptr)
                insertWithoutLookup(expanded, entry);
        }
        slots.swap(expanded);
    }
};

#endif
