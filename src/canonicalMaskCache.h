#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "compilerAttributes.h"

/** Estimate container storage, excluding allocator headers and padding. */
template<typename Map>
[[nodiscard]] std::uint64_t unorderedRetainedBytes(const Map &map) noexcept
{
    // A default empty map commonly owns no bucket allocation.
    const std::uint64_t buckets = map.bucket_count() > 1
        ? map.bucket_count() * sizeof(void *) : 0;
    return buckets + map.size() * (
        sizeof(typename Map::value_type) + 2 * sizeof(void *)
    );
}

/**
 * Exact mask-to-ID lookup with an optional two-word flat front table.
 *
 * Each generation starts in the unordered map; small local working sets keep
 * its cheaper hash. On the 512th distinct mask, eligible domains promote all
 * entries into a flat table. Its allocation is capped at 65,536 slots.
 * Overflow retains the original unordered map semantics: DAG finalization
 * requires every enumerated mask to remain present. Domains wider than 128
 * bits use only the unordered map. Exact canonicalization is unchanged.
 */
template<typename Mask>
class canonicalMaskCache
{
public:
    enum class policy { unordered, flat };
    using mapped_type = std::pair<int, int>;
    struct value_type { mapped_type second{-1, 0}; };
    using const_iterator = const value_type *;
    static constexpr std::size_t flatPromotionThreshold = 512;
    static constexpr std::size_t maximumFlatCapacity = 65536;

    explicit canonicalMaskCache(policy selected = policyFromEnvironment()):
        selectedPolicy(selected) {}

    [[nodiscard]] static policy policyFromEnvironment() noexcept
    {
        const char *setting = std::getenv("ASSEMBLY_CANONICAL_MASK_CACHE");
        return setting != nullptr && std::string_view(setting) == "unordered"
            ? policy::unordered : policy::flat;
    }

    /** Clear before reconfiguring Mask's domain, as for a mask-keyed map. */
    void clear() noexcept
    {
        overflow.clear();
        for (slot &entry : slots) entry.second.first = -1;
        flatSize = 0;
        flatActive = false;
    }

    void setPolicy(policy selected)
    {
        if (!empty()) throw std::logic_error("mask cache policy requires an empty cache");
        selectedPolicy = selected;
    }

    [[nodiscard]] bool usesFlatStorage() const noexcept
    {
        return flatActive;
    }

    [[nodiscard]] std::size_t size() const noexcept
    {
        return flatSize + overflow.size();
    }

    [[nodiscard]] bool empty() const noexcept { return size() == 0; }
    [[nodiscard]] const_iterator end() const noexcept { return nullptr; }
    [[nodiscard]] std::size_t flatCapacity() const noexcept { return slots.size(); }
    [[nodiscard]] std::size_t overflowSize() const noexcept { return overflow.size(); }

    [[nodiscard]] const_iterator find(const Mask &mask) const noexcept
    {
        if (flatActive)
        {
            const auto [low, high] = words(mask);
            const slot &entry = slots[findSlot(slots, low, high)];
            if (entry.second.first >= 0) return &entry;
            if (overflow.empty()) return end();
        }
        const auto found = overflow.find(mask);
        return found == overflow.end() ? end() : &found->second;
    }

    void emplace(const Mask &mask, mapped_type canonical)
    {
        if (canonical.first < 0)
            throw std::invalid_argument("canonical mask ID must be nonnegative");
        if (flatActive)
        {
            const auto [low, high] = words(mask);
            std::size_t position = findSlot(slots, low, high);
            if (slots[position].second.first >= 0) return;
            if (flatSize < maximumEntries(slots.size()))
            {
                insert(position, low, high, canonical);
                return;
            }
            if (slots.size() < maximumFlatCapacity)
            {
                grow(slots.size() * 2);
                position = findSlot(slots, low, high);
                insert(position, low, high, canonical);
                return;
            }
        }
        else if (
            selectedPolicy == policy::flat && Mask::activeWordCount() <= 2 &&
            overflow.size() >= flatPromotionThreshold - 1
        )
        {
            if (overflow.find(mask) != overflow.end()) return;
            promote(mask, canonical);
            return;
        }
        overflow.emplace(mask, value_type{canonical});
    }

    /**
     * Exact flat capacity plus estimated unordered nodes/buckets and the wide
     * words retained by its distinct mask keys. Excludes arena slack, allocator
     * headers and shared references to those words elsewhere in the search.
     */
    [[nodiscard]] std::uint64_t retainedBytes() const noexcept
    {
        std::uint64_t bytes = slots.capacity() * sizeof(slot) +
            unorderedRetainedBytes(overflow);
        if (Mask::activeWordCount() > 1)
        {
            for (const auto &entry : overflow)
            {
                if (entry.first != 0)
                    bytes += sizeof(std::size_t) +
                        Mask::activeWordCount() * sizeof(std::uint64_t);
            }
        }
        return bytes;
    }

private:
    struct slot : value_type
    {
        std::uint64_t low = 0;
        std::uint64_t high = 0;
    };
    static_assert(sizeof(slot) == 24);
    static constexpr std::size_t initialFlatCapacity = 1024;
    policy selectedPolicy;
    bool flatActive = false;
    std::vector<slot> slots;
    std::unordered_map<Mask, value_type> overflow;
    std::size_t flatSize = 0;

    [[nodiscard]] static std::pair<std::uint64_t, std::uint64_t> words(
        const Mask &mask
    ) noexcept
    {
        return {
            Mask::activeWordCount() == 0 ? 0 : mask.activeWord(0),
            Mask::activeWordCount() < 2 ? 0 : mask.activeWord(1)
        };
    }

    [[nodiscard]] static std::size_t hashWords(
        std::uint64_t low, std::uint64_t high
    ) noexcept
    {
        std::uint64_t hash = low ^ (high + 0x9e3779b97f4a7c15ULL +
            (low << 6) + (low >> 2));
        hash ^= hash >> 30;
        hash *= 0xbf58476d1ce4e5b9ULL;
        hash ^= hash >> 27;
        hash *= 0x94d049bb133111ebULL;
        return hash ^ (hash >> 31);
    }

    [[nodiscard]] static std::size_t maximumEntries(std::size_t capacity) noexcept
    {
        return capacity * 7 / 10;
    }

    [[nodiscard]] static std::size_t findSlot(
        const std::vector<slot> &table,
        std::uint64_t low,
        std::uint64_t high
    ) noexcept
    {
        const std::size_t indexMask = table.size() - 1;
        std::size_t index = hashWords(low, high) & indexMask;
        while (table[index].second.first >= 0 &&
            (table[index].low != low || table[index].high != high))
        {
            index = (index + 1) & indexMask;
        }
        return index;
    }

    void insert(
        std::size_t position, std::uint64_t low, std::uint64_t high,
        mapped_type canonical
    ) noexcept
    {
        slots[position] = slot{{canonical}, low, high};
        ++flatSize;
    }

    /**
     * Allocate before mutation, then copy only scalar words and IDs. A failed
     * allocation leaves every existing entry and the inactive state intact;
     * the incoming mask is inserted only after allocation succeeds. Clearing
     * a generation retains empty slots which can be reused without allocating.
     */
    PARALLELASSEMBLYCPP_NOINLINE void promote(
        const Mask &mask, mapped_type canonical
    )
    {
        decltype(overflow) previousGeneration;
        if (slots.empty()) grow(initialFlatCapacity);
        for (const auto &entry : overflow)
        {
            const auto [low, high] = words(entry.first);
            insert(findSlot(slots, low, high), low, high, entry.second.second);
        }
        const auto [low, high] = words(mask);
        insert(findSlot(slots, low, high), low, high, canonical);
        flatActive = true;
        // Destroy old nodes and bucket storage, releasing owning wide masks.
        overflow.swap(previousGeneration);
    }

    PARALLELASSEMBLYCPP_NOINLINE void grow(std::size_t capacity)
    {
        std::vector<slot> expanded(capacity);
        for (const slot &entry : slots)
        {
            if (entry.second.first >= 0)
                expanded[findSlot(expanded, entry.low, entry.high)] = entry;
        }
        slots.swap(expanded);
    }
};
