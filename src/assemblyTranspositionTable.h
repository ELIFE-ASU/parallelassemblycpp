#ifndef ASSEMBLY_TRANSPOSITION_TABLE_H
#define ASSEMBLY_TRANSPOSITION_TABLE_H

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <cstdlib>
#include <string_view>
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    #include <chrono>
#endif
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

#include "compilerAttributes.h"

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

    PARALLELASSEMBLYCPP_NOINLINE result insertMiss(
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
 * Exact process-shared L2 with independent lock shards. A contiguous hash
 * index rejects nonmatches without loading the parallel key-pointer index.
 * Scores live with exact keys; growth reads hashes without chasing pointers.
 * Admission is optional: rejected misses never imply domination. Existing
 * entries remain available for exact lookup and monotonic score improvement.
 */
class sharedAssemblyTranspositionTable
{
public:
    static constexpr std::size_t shardCount = 64;
    enum class policy { shared, local, selective };

    /** Experiment controls, read once before workers enter the search. */
    static policy policyFromEnvironment()
    {
#if defined(_MSC_VER)
        char *buffer = nullptr;
        std::size_t length = 0;
        const auto status = _dupenv_s(
            &buffer, &length, "PARALLELASSEMBLYCPP_SHARED_CACHE_POLICY"
        );
        const std::unique_ptr<char, decltype(&std::free)> owner(buffer, &std::free);
        if (status != 0)
            throw std::runtime_error("could not read shared cache policy");
        const char *value = owner.get();
#else
        const char *value = std::getenv("PARALLELASSEMBLYCPP_SHARED_CACHE_POLICY");
#endif
        if (value == nullptr || std::string_view(value) == "shared")
            return policy::shared;
        if (std::string_view(value) == "local") return policy::local;
        if (std::string_view(value) == "selective") return policy::selective;
        throw std::invalid_argument("invalid shared cache policy");
    }

    static std::size_t maxBytesFromEnvironment()
    {
#if defined(_MSC_VER)
        char *buffer = nullptr;
        std::size_t length = 0;
        const auto status = _dupenv_s(
            &buffer, &length, "PARALLELASSEMBLYCPP_SHARED_CACHE_BYTES"
        );
        const std::unique_ptr<char, decltype(&std::free)> owner(buffer, &std::free);
        if (status != 0)
            throw std::runtime_error("could not read shared cache byte budget");
        const char *value = owner.get();
#else
        const char *value = std::getenv("PARALLELASSEMBLYCPP_SHARED_CACHE_BYTES");
#endif
        if (value == nullptr) return 0;
        const std::string_view text(value);
        std::size_t bytes = 0;
        const auto parsed = std::from_chars(
            text.data(), text.data() + text.size(), bytes
        );
        if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
            throw std::invalid_argument("invalid shared cache byte budget");
        return bytes;
    }

    /**
     * Worker indices own unsynchronised arenas; concurrent callers must use
     * distinct indices. The upstream must outlive the table and support
     * concurrent refills. maxBytes caps retained entry/key bytes, split evenly
     * over shards; zero means unlimited. Index arrays and arena slack are
     * reported separately and are not included in this admission budget.
     */
    explicit sharedAssemblyTranspositionTable(
        std::size_t workerCount = 0,
        std::pmr::memory_resource *workerUpstream =
            std::pmr::new_delete_resource(),
        policy admissionPolicy = policy::shared,
        std::size_t maxBytes = 0
    ):
        allocationResource(workerUpstream),
        admissionPolicyValue(admissionPolicy),
        bounded(maxBytes != 0),
        shardByteBudget(maxBytes / shardCount)
    {
        workerPools.reserve(workerCount);
        for (std::size_t worker = 0; worker < workerCount; ++worker)
            workerPools.push_back(std::make_unique<
                std::pmr::monotonic_buffer_resource
            >(64 * 1024, &allocationResource));
    }

    ~sharedAssemblyTranspositionTable()
    {
        for (shard &selected : shards)
        {
            // Production workers publish only arena-owned, trivially
            // destructible keys. Release those arenas in bulk without a
            // pointer-chasing pass over every occupied slot at shutdown.
            if (selected.directAllocatedBytes == 0) continue;
            for (entry *stored : selected.activeKeys())
            {
                if (stored != nullptr && stored->pooled == 0)
                {
                    const std::size_t bytes = entryBytes(stored->length);
                    std::pmr::new_delete_resource()->deallocate(
                        stored, bytes, alignof(entry)
                    );
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

    /**
     * Completed lookups satisfy hits + misses == lock acquisitions,
     * hits == pruned + updated, and misses == admissions + rejections.
     * allocatedBytes is retained entry/key storage; slotBytes includes inline
     * and expanded indices, arenaAllocatedBytes includes upstream arena slack
     * and directly allocated keys. Timings are collected only in telemetry
     * builds. Growth duration includes allocation, rehash and old-index free
     * under the lock; maxGrowthNanoseconds is the longest individual growth.
     */
    struct statistics
    {
        std::uint64_t hitCount = 0;
        std::uint64_t missCount = 0;
        std::uint64_t collisionChainSteps = 0;
        std::uint64_t allocatedBytes = 0;
        std::uint64_t lockAcquisitionCount = 0;
        std::uint64_t lockWaitCount = 0;
        std::uint64_t lockWaitNanoseconds = 0;
        std::uint64_t admissionCount = 0;
        std::uint64_t admissionRejectionCount = 0;
        std::uint64_t prunedHitCount = 0;
        std::uint64_t updatedHitCount = 0;
        std::uint64_t slotBytes = 0;
        std::uint64_t arenaAllocatedBytes = 0;
        std::uint64_t admissionFilterBytes = 0;
        std::uint64_t growthCount = 0;
        std::uint64_t rehashedEntries = 0;
        std::uint64_t growthNanoseconds = 0;
        std::uint64_t maxGrowthNanoseconds = 0;
    };

    [[nodiscard]] bool lookupEnabled() const noexcept
    {
        return admissionPolicyValue != policy::local;
    }

    assemblyTranspositionTable::result consider(
        std::span<const int> key, int sumDupBonds
    )
    {
        return considerWithBest(key, sumDupBonds).outcome;
    }

    consideration considerWithBest(std::span<const int> key, int sumDupBonds)
    {
        return considerWithResource(
            key, sumDupBonds, *std::pmr::new_delete_resource(), false
        );
    }

    consideration considerWithBestForWorker(
        std::span<const int> key, int sumDupBonds, std::size_t workerIndex
    )
    {
        if (workerIndex >= workerPools.size())
            throw std::out_of_range("shared-table worker index is invalid");
        return considerWithResource(
            key, sumDupBonds, *workerPools[workerIndex], true
        );
    }

    [[nodiscard]] std::size_t workerCount() const noexcept
    {
        return workerPools.size();
    }

    [[nodiscard]] std::size_t size() const
    {
        std::size_t result = 0;
        for (const shard &selected : shards)
        {
            std::lock_guard lock(selected.mutex);
            result += selected.sizeValue;
        }
        return result;
    }

    [[nodiscard]] statistics stats() const
    {
        statistics result;
        for (const shard &selected : shards)
        {
            std::lock_guard lock(selected.mutex);
            const statistics &source = selected.counters;
            result.hitCount += source.hitCount;
            result.missCount += source.missCount;
            result.collisionChainSteps += source.collisionChainSteps;
            result.allocatedBytes += source.allocatedBytes;
            result.lockAcquisitionCount += source.lockAcquisitionCount;
            result.lockWaitCount += source.lockWaitCount;
            result.lockWaitNanoseconds += source.lockWaitNanoseconds;
            result.admissionCount += source.admissionCount;
            result.admissionRejectionCount += source.admissionRejectionCount;
            result.prunedHitCount += source.prunedHitCount;
            result.updatedHitCount += source.updatedHitCount;
            result.growthCount += source.growthCount;
            result.rehashedEntries += source.rehashedEntries;
            result.growthNanoseconds += source.growthNanoseconds;
            result.maxGrowthNanoseconds = std::max(
                result.maxGrowthNanoseconds, source.maxGrowthNanoseconds
            );
            result.slotBytes += sizeof(selected.initialKeys) +
                sizeof(selected.initialHashes) +
                selected.expandedKeys.capacity() * sizeof(entry *) +
                selected.expandedHashes.capacity() * sizeof(std::uint32_t);
            result.admissionFilterBytes +=
                selected.admissionFilter.capacity() * sizeof(std::uint32_t);
            result.arenaAllocatedBytes += selected.directAllocatedBytes;
        }
        result.arenaAllocatedBytes += allocationResource.bytes.load(
            std::memory_order_relaxed
        );
        return result;
    }

private:
    static constexpr std::size_t minimumShardCapacity = 1024;
    static constexpr std::size_t admissionFilterCapacity = 4096;

    struct alignas(int) entry
    {
        std::uint32_t length;
        int bestSumDupBonds;
        std::uint32_t pooled;
    };
    static_assert(sizeof(entry) == 12);
    static_assert(std::is_trivially_destructible_v<entry>);

    struct entryDeleter
    {
        std::pmr::memory_resource *resource = nullptr;
        std::size_t bytes = 0;
        void operator()(entry *value) const noexcept
        {
            if (value != nullptr)
                resource->deallocate(value, bytes, alignof(entry));
        }
    };
    using ownedEntry = std::unique_ptr<entry, entryDeleter>;

    struct entryFindResult
    {
        entry *value = nullptr;
        std::size_t index = 0;
        std::uint64_t collisionChainSteps = 0;
    };

    struct alignas(64) shard
    {
        mutable std::mutex mutex;
        std::array<std::uint32_t, minimumShardCapacity> initialHashes{};
        std::array<entry *, minimumShardCapacity> initialKeys{};
        std::vector<std::uint32_t> expandedHashes;
        std::vector<entry *> expandedKeys;
        // A bounded, direct-mapped first-sighting filter, allocated only for
        // selective mode. Fingerprints control admission, never equality.
        std::vector<std::uint32_t> admissionFilter;
        std::size_t sizeValue = 0;
        std::uint64_t directAllocatedBytes = 0;
        statistics counters;

        std::span<std::uint32_t> activeHashes() noexcept
        {
            if (expandedHashes.empty()) return initialHashes;
            return expandedHashes;
        }

        std::span<entry *> activeKeys() noexcept
        {
            if (expandedKeys.empty()) return initialKeys;
            return expandedKeys;
        }
    };

    class measuredResource final : public std::pmr::memory_resource
    {
    public:
        std::atomic<std::uint64_t> bytes{0};
        explicit measuredResource(std::pmr::memory_resource *resource):
            upstream(resource)
        {
            if (upstream == nullptr)
                throw std::invalid_argument("shared-table worker upstream is null");
        }
    private:
        std::pmr::memory_resource *upstream;
        void *do_allocate(std::size_t size, std::size_t alignment) override
        {
            void *result = upstream->allocate(size, alignment);
            bytes.fetch_add(size, std::memory_order_relaxed);
            return result;
        }
        void do_deallocate(
            void *memory, std::size_t size, std::size_t alignment
        ) override
        {
            upstream->deallocate(memory, size, alignment);
            bytes.fetch_sub(size, std::memory_order_relaxed);
        }
        bool do_is_equal(const std::pmr::memory_resource &other) const
            noexcept override
        {
            return this == &other;
        }
    };

    measuredResource allocationResource;
    std::vector<std::unique_ptr<std::pmr::monotonic_buffer_resource>> workerPools;
    std::array<shard, shardCount> shards;
    policy admissionPolicyValue;
    bool bounded;
    std::size_t shardByteBudget;

    consideration considerWithResource(
        std::span<const int> key,
        int sumDupBonds,
        std::pmr::memory_resource &resource,
        bool pooled
    )
    {
        if (key.size() > std::numeric_limits<std::uint32_t>::max() ||
            key.size_bytes() >
                std::numeric_limits<std::size_t>::max() - sizeof(entry))
            throw std::length_error("assembly-state key is too long");
        if (!lookupEnabled())
            return {assemblyTranspositionTable::result::inserted, sumDupBonds};

        const std::uint32_t hash = assemblyTranspositionTable::keyHash(key);
        shard &selected = shards[shardIndex(hash)];
        bool lockWaited = false;
        std::uint64_t lockWaitNanoseconds = 0;
        std::unique_lock lock = lockShard(
            selected, lockWaited, lockWaitNanoseconds
        );
        const entryFindResult found = find(selected, key, hash);
        selected.counters.collisionChainSteps += found.collisionChainSteps;
        if (found.value != nullptr)
        {
            ++selected.counters.hitCount;
            recordCompletedLock(selected, lockWaited, lockWaitNanoseconds);
            if (sumDupBonds <= found.value->bestSumDupBonds)
            {
                ++selected.counters.prunedHitCount;
                return {
                    assemblyTranspositionTable::result::dominated,
                    found.value->bestSumDupBonds
                };
            }
            ++selected.counters.updatedHitCount;
            found.value->bestSumDupBonds = sumDupBonds;
            return {assemblyTranspositionTable::result::improved, sumDupBonds};
        }

        const std::size_t bytes = entryBytes(key.size());
        const bool withinBudget = !bounded ||
            bytes <= shardByteBudget - selected.counters.allocatedBytes;
        if (!withinBudget || !admit(selected, hash))
        {
            ++selected.counters.missCount;
            ++selected.counters.admissionRejectionCount;
            recordCompletedLock(selected, lockWaited, lockWaitNanoseconds);
            return {assemblyTranspositionTable::result::inserted, sumDupBonds};
        }

        std::size_t insertionIndex = found.index;
        if (selected.sizeValue >= maximumShardEntries(selected.activeHashes().size()))
        {
            grow(selected);
            insertionIndex = emptySlot(selected.activeHashes(), hash);
        }
        // Allocate only after absence and admission are established. Neither
        // rejected misses nor racing hits consume monotonic arena space.
        ownedEntry prepared = prepareEntry(key, sumDupBonds, resource, pooled);
        selected.activeKeys()[insertionIndex] = prepared.release();
        // The low bits already selected the shard; setting bit zero cannot
        // change the slot index. Zero remains an unambiguous empty marker.
        selected.activeHashes()[insertionIndex] = hash | 1U;
        ++selected.sizeValue;
        ++selected.counters.missCount;
        ++selected.counters.admissionCount;
        selected.counters.allocatedBytes += bytes;
        if (!pooled) selected.directAllocatedBytes += bytes;
        recordCompletedLock(selected, lockWaited, lockWaitNanoseconds);
        return {assemblyTranspositionTable::result::inserted, sumDupBonds};
    }

    bool admit(shard &selected, std::uint32_t hash)
    {
        if (admissionPolicyValue != policy::selective) return true;
        if (selected.admissionFilter.empty())
            selected.admissionFilter.resize(admissionFilterCapacity);
        std::uint32_t &seen = selected.admissionFilter[
            slotIndex(hash, admissionFilterCapacity - 1)
        ];
        const std::uint32_t fingerprint = hash | 1U;
        if (seen == fingerprint) return true;
        seen = fingerprint;
        return false;
    }

    static bool keysEqual(const entry &stored, std::span<const int> key) noexcept
    {
        return stored.length == key.size() && std::equal(
            key.begin(), key.end(),
            reinterpret_cast<const int *>(std::addressof(stored) + 1)
        );
    }

    static entryFindResult find(
        shard &selected, std::span<const int> key, std::uint32_t hash
    ) noexcept
    {
        entryFindResult result;
        const std::span<std::uint32_t> hashes = selected.activeHashes();
        const std::size_t mask = hashes.size() - 1;
        std::size_t index = slotIndex(hash, mask);
        const std::uint32_t storedHash = hash | 1U;
        while (hashes[index] != 0)
        {
            if (hashes[index] == storedHash)
            {
                entry *candidate = selected.activeKeys()[index];
                if (keysEqual(*candidate, key))
                {
                    result.value = candidate;
                    result.index = index;
                    return result;
                }
            }
            ++result.collisionChainSteps;
            index = (index + 1) & mask;
        }
        result.index = index;
        return result;
    }

    static std::size_t maximumShardEntries(std::size_t capacity) noexcept
    {
        return capacity - (capacity + 4) / 5;
    }

    static std::size_t slotIndex(std::uint32_t hash, std::size_t mask) noexcept
    {
        return (static_cast<std::size_t>(hash) / shardCount) & mask;
    }

    static std::size_t emptySlot(
        std::span<const std::uint32_t> destination, std::uint32_t hash
    ) noexcept
    {
        const std::size_t mask = destination.size() - 1;
        std::size_t index = slotIndex(hash, mask);
        while (destination[index] != 0)
            index = (index + 1) & mask;
        return index;
    }

    static void grow(shard &selected)
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        const auto started = std::chrono::steady_clock::now();
#endif
        {
            const std::span<std::uint32_t> hashes = selected.activeHashes();
            const std::span<entry *> keys = selected.activeKeys();
            const std::size_t maximumCapacity = std::min(
                selected.expandedHashes.max_size(),
                selected.expandedKeys.max_size()
            );
            if (hashes.size() > maximumCapacity / 4)
                throw std::length_error("shared transposition shard is too large");
            // Allocate both indices before mutating the shard. If either
            // allocation fails, every existing key remains reachable.
            std::vector<std::uint32_t> expandedHashes(hashes.size() * 4);
            std::vector<entry *> expandedKeys(hashes.size() * 4);
            for (std::size_t index = 0; index < hashes.size(); ++index)
            {
                if (hashes[index] == 0) continue;
                const std::size_t destination =
                    emptySlot(expandedHashes, hashes[index]);
                expandedHashes[destination] = hashes[index];
                expandedKeys[destination] = keys[index];
            }
            selected.expandedHashes.swap(expandedHashes);
            selected.expandedKeys.swap(expandedKeys);
        }
        ++selected.counters.growthCount;
        selected.counters.rehashedEntries += selected.sizeValue;
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        const auto elapsed = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started
            ).count()
        );
        selected.counters.growthNanoseconds += elapsed;
        selected.counters.maxGrowthNanoseconds = std::max(
            selected.counters.maxGrowthNanoseconds, elapsed
        );
#endif
    }

    static std::unique_lock<std::mutex> lockShard(
        shard &selected, bool &waited, std::uint64_t &waitNanoseconds
    )
    {
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        std::unique_lock<std::mutex> lock(selected.mutex, std::try_to_lock);
        if (lock.owns_lock()) return lock;
        const auto started = std::chrono::steady_clock::now();
        lock.lock();
        waited = true;
        waitNanoseconds = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - started
            ).count()
        );
        return lock;
#else
        static_cast<void>(waited);
        static_cast<void>(waitNanoseconds);
        return std::unique_lock<std::mutex>(selected.mutex);
#endif
    }

    static void recordCompletedLock(
        shard &selected, bool waited, std::uint64_t waitNanoseconds
    ) noexcept
    {
        ++selected.counters.lockAcquisitionCount;
        if (waited)
        {
            ++selected.counters.lockWaitCount;
            selected.counters.lockWaitNanoseconds += waitNanoseconds;
        }
    }

    static ownedEntry prepareEntry(
        std::span<const int> key, int score,
        std::pmr::memory_resource &resource, bool pooled
    )
    {
        const std::size_t bytes = entryBytes(key.size());
        void *memory = resource.allocate(bytes, alignof(entry));
        ownedEntry prepared(
            std::construct_at(static_cast<entry *>(memory), entry{
                static_cast<std::uint32_t>(key.size()),
                score,
                static_cast<std::uint32_t>(pooled)
            }),
            entryDeleter{std::addressof(resource), bytes}
        );
        std::uninitialized_copy(
            key.begin(), key.end(), reinterpret_cast<int *>(prepared.get() + 1)
        );
        return prepared;
    }

    static std::size_t entryBytes(std::size_t keyLength) noexcept
    {
        return sizeof(entry) + sizeof(int) * keyLength;
    }

    static std::size_t shardIndex(std::uint32_t hash) noexcept
    {
        static_assert((shardCount & (shardCount - 1)) == 0);
        return static_cast<std::size_t>(hash) & (shardCount - 1);
    }
};

#endif
