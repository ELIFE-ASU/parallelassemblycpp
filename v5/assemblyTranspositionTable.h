#ifndef ASSEMBLY_TRANSPOSITION_TABLE_H
#define ASSEMBLY_TRANSPOSITION_TABLE_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <memory_resource>
#include <mutex>
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

    /** Hash a borrowed key before entering an optional external lock. */
    [[nodiscard]] static std::uint32_t keyHash(
        std::span<const int> key
    ) noexcept
    {
        return hashKey(key);
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

/**
 * @brief Exact process-shared L2 table split across independent lock shards.
 *
 * A worker computes and canonicalises its complete key before calling this
 * table. Hashing, key copying, and node allocation happen before the shard
 * lock; the fixed-bucket critical section is limited to one exact checked
 * pointer link, score improvement, or dominated lookup and never resizes.
 */
class sharedAssemblyTranspositionTable
{
public:
    static constexpr std::size_t shardCount = 64;

    explicit sharedAssemblyTranspositionTable(std::size_t workerCount = 0)
    {
        workerPools.reserve(workerCount);
        for (std::size_t worker = 0; worker < workerCount; ++worker)
        {
            workerPools.push_back(std::make_unique<
                std::pmr::monotonic_buffer_resource
            >(64 * 1024, std::pmr::new_delete_resource()));
        }
    }

    ~sharedAssemblyTranspositionTable()
    {
        for (shard &selected : shards)
        {
            for (entry *head : selected.buckets)
            {
                while (head != nullptr)
                {
                    entry *next = head->next;
                    const bool pooled = head->pooled != 0;
                    const std::size_t bytes =
                        sizeof(entry) + sizeof(int) * head->length;
                    std::destroy_at(head);
                    if (!pooled)
                    {
                        std::pmr::new_delete_resource()->deallocate(
                            head,
                            bytes,
                            alignof(entry)
                        );
                    }
                    head = next;
                }
            }
        }
    }

    sharedAssemblyTranspositionTable(
        const sharedAssemblyTranspositionTable &
    ) = delete;
    sharedAssemblyTranspositionTable &operator=(
        const sharedAssemblyTranspositionTable &
    ) = delete;

    struct consideration
    {
        assemblyTranspositionTable::result outcome;
        int bestSumDupBonds;
    };

    assemblyTranspositionTable::result consider(
        std::span<const int> key,
        int sumDupBonds
    )
    {
        return considerWithBest(key, sumDupBonds).outcome;
    }

    consideration considerWithBest(
        std::span<const int> key,
        int sumDupBonds
    )
    {
        return considerWithResource(
            key,
            sumDupBonds,
            *std::pmr::new_delete_resource(),
            false
        );
    }

    consideration considerWithBestForWorker(
        std::span<const int> key,
        int sumDupBonds,
        std::size_t workerIndex
    )
    {
        // A worker index selects an unsynchronised monotonic arena. Concurrent
        // callers must therefore use distinct indices; the OpenMP wiring
        // assigns its stable thread index for the lifetime of the search.
        if (workerIndex >= workerPools.size())
            throw std::out_of_range("shared-table worker index is invalid");
        return considerWithResource(
            key,
            sumDupBonds,
            *workerPools[workerIndex],
            true
        );
    }

    [[nodiscard]] std::size_t workerCount() const noexcept
    {
        return workerPools.size();
    }

    [[nodiscard]] std::size_t size() const
    {
        std::size_t result = 0;
        for (std::size_t index = 0; index < shardCount; ++index)
        {
            const shard &selected = shards[index];
            std::lock_guard lock(selected.mutex);
            result += selected.sizeValue;
        }
        return result;
    }

private:
    static constexpr std::size_t bucketCount = 1024;

    struct alignas(int) entry
    {
        entry *next;
        std::uint32_t hash;
        std::uint32_t length;
        int bestSumDupBonds;
        std::uint32_t pooled;
    };

    static_assert(sizeof(entry) % alignof(int) == 0);

    struct entryDeleter
    {
        std::pmr::memory_resource *resource = nullptr;
        std::size_t bytes = 0;

        void operator()(entry *value) const noexcept
        {
            if (value == nullptr) return;
            std::destroy_at(value);
            resource->deallocate(value, bytes, alignof(entry));
        }
    };

    using ownedEntry = std::unique_ptr<entry, entryDeleter>;

    struct alignas(64) shard
    {
        mutable std::mutex mutex;
        std::array<entry *, bucketCount> buckets{};
        std::size_t sizeValue = 0;
    };

    // Pools are declared before shards so their storage remains alive until
    // after every published entry header has been destroyed.
    std::vector<
        std::unique_ptr<std::pmr::monotonic_buffer_resource>
    > workerPools;
    std::array<shard, shardCount> shards;

    consideration considerWithResource(
        std::span<const int> key,
        int sumDupBonds,
        std::pmr::memory_resource &resource,
        bool pooled
    )
    {
        if (key.size() > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("assembly-state key is too long");

        const std::uint32_t hash = assemblyTranspositionTable::keyHash(key);
        // The worker arena makes preparation a bump allocation. Preparing
        // before the lock avoids a second mutex round trip on every L2 miss;
        // a cross-worker hit simply drops the unused prepared node afterward.
        ownedEntry prepared = prepareEntry(
            key,
            sumDupBonds,
            hash,
            resource,
            pooled
        );
        shard &selected = shards[shardIndex(hash)];
        {
            std::lock_guard lock(selected.mutex);
            if (entry *existing = find(selected, key, hash))
                return updateExisting(*existing, sumDupBonds);
            const std::size_t bucket = bucketIndex(hash);
            prepared->next = selected.buckets[bucket];
            selected.buckets[bucket] = prepared.release();
            ++selected.sizeValue;
        }
        return {
            assemblyTranspositionTable::result::inserted,
            sumDupBonds
        };
    }

    static const int *keyValues(const entry &value) noexcept
    {
        return reinterpret_cast<const int *>(std::addressof(value) + 1);
    }

    static bool keysEqual(
        const entry &stored,
        std::span<const int> candidate
    ) noexcept
    {
        return stored.length == candidate.size() && std::equal(
            candidate.begin(),
            candidate.end(),
            keyValues(stored)
        );
    }

    static entry *find(
        shard &selected,
        std::span<const int> key,
        std::uint32_t hash
    ) noexcept
    {
        entry *candidate = selected.buckets[bucketIndex(hash)];
        while (candidate != nullptr)
        {
            if (
                candidate->hash == hash && keysEqual(*candidate, key)
            ) return candidate;
            candidate = candidate->next;
        }
        return nullptr;
    }

    static consideration updateExisting(entry &existing, int score) noexcept
    {
        if (score <= existing.bestSumDupBonds)
        {
            return {
                assemblyTranspositionTable::result::dominated,
                existing.bestSumDupBonds
            };
        }
        existing.bestSumDupBonds = score;
        return {assemblyTranspositionTable::result::improved, score};
    }

    static ownedEntry prepareEntry(
        std::span<const int> key,
        int score,
        std::uint32_t hash,
        std::pmr::memory_resource &resource,
        bool pooled
    )
    {
        if (
            key.size_bytes() >
            std::numeric_limits<std::size_t>::max() - sizeof(entry)
        ) throw std::length_error("assembly-state key allocation overflow");
        const std::size_t bytes = sizeof(entry) + key.size_bytes();
        void *memory = resource.allocate(bytes, alignof(entry));
        ownedEntry prepared(
            std::construct_at(
                static_cast<entry *>(memory),
                entry{
                    nullptr,
                    hash,
                    static_cast<std::uint32_t>(key.size()),
                    score,
                    static_cast<std::uint32_t>(pooled)
                }
            ),
            entryDeleter{std::addressof(resource), bytes}
        );
        std::uninitialized_copy(
            key.begin(),
            key.end(),
            reinterpret_cast<int *>(prepared.get() + 1)
        );
        return prepared;
    }

    static std::size_t shardIndex(std::uint32_t hash) noexcept
    {
        static_assert((shardCount & (shardCount - 1)) == 0);
        return static_cast<std::size_t>(hash) & (shardCount - 1);
    }

    static std::size_t bucketIndex(std::uint32_t hash) noexcept
    {
        static_assert((bucketCount & (bucketCount - 1)) == 0);
        return (
            static_cast<std::size_t>(hash) /
            shardCount
        ) & (bucketCount - 1);
    }
};

#endif
