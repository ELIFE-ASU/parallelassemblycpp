#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../src/activeWordMask.h"
#include "../src/canonicalMaskCache.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

using Cache = canonicalMaskCache<EdgeMask>;

void testWidth(std::size_t width, Cache::policy policy)
{
    EdgeMask::configure(width);
    Cache cache(policy);
    std::unordered_map<EdgeMask, std::pair<int, int>> expected;
    assert(cache.empty());
    assert(cache.retainedBytes() == 0);
    assert(!cache.usesFlatStorage());
    {
        EdgeMask mask;
        cache.emplace(mask, {0, 2});
        expected.emplace(mask, std::pair{0, 2});
        const auto initial = cache.find(mask);
        assert(initial != cache.end());
        assert(initial->second.first == 0);
        for (std::size_t bit = 0; bit < width; ++bit)
        {
            mask.set(bit);
            const std::pair<int, int> value{static_cast<int>(bit + 1), 3};
            cache.emplace(mask, value);
            expected.emplace(mask, value);
        }
        assert(cache.size() == expected.size());
        for (const auto &[key, value] : expected)
        {
            const auto found = cache.find(key);
            assert(found != cache.end());
            assert(found->second == value);
            cache.emplace(key, {999, 999});
            const auto repeated = cache.find(key);
            assert(repeated != cache.end());
            assert(repeated->second == value);
        }
        if (width > 1)
        {
            EdgeMask absent;
            absent.set(width - 1);
            assert(cache.find(absent) == cache.end());
        }
        bool rejected = false;
        try { cache.setPolicy(Cache::policy::unordered); }
        catch (const std::logic_error &) { rejected = true; }
        assert(rejected);
    }
    assert(cache.retainedBytes() > 0);
    expected.clear();
    const auto capacity = cache.flatCapacity();
    cache.clear();
    assert(cache.empty());
    assert(cache.flatCapacity() == capacity);
    assert(cache.retainedBytes() >= capacity * 24);
    EdgeMask empty;
    assert(cache.find(empty) == cache.end());
    cache.emplace(empty, {42, 1});
    const auto reinserted = cache.find(empty);
    assert(reinserted != cache.end());
    assert(reinserted->second.first == 42);
    cache.clear();
}

EdgeMask numberedMask(std::uint64_t value)
{
    const std::size_t wordCount = EdgeMask::activeWordCount();
    if (wordCount == 0 || wordCount > 2)
        throw std::logic_error("numbered test mask requires one or two words");
    std::array<std::uint64_t, 2> words{value, 0};
    if (wordCount == 2)
    {
        words[0] = value & 255;
        words[1] = value >> 8;
    }
    return EdgeMask::fromActiveWords(words.data());
}

void requireNumberedEntry(const Cache &cache, std::uint64_t value)
{
    const EdgeMask mask = numberedMask(value);
    const auto found = cache.find(mask);
    assert(found != cache.end());
    assert(found->second.first == static_cast<int>(value));
    assert(found->second.second == 7);
}

void testPromotionAndGenerationReset(std::size_t width)
{
    EdgeMask::configure(width);
    Cache cache(Cache::policy::flat);
    for (std::uint64_t i = 0; i < Cache::flatPromotionThreshold - 1; ++i)
    {
        const EdgeMask mask = numberedMask(i);
        cache.emplace(mask, {static_cast<int>(i), 7});
    }
    assert(cache.size() == Cache::flatPromotionThreshold - 1);
    assert(!cache.usesFlatStorage());
    assert(cache.flatCapacity() == 0);
    // Duplicate insertions at the threshold preserve values and stay in the
    // small-set backend; only a new distinct mask can trigger promotion.
    for (std::uint64_t i = 0; i < Cache::flatPromotionThreshold - 1; ++i)
    {
        const EdgeMask mask = numberedMask(i);
        cache.emplace(mask, {9999, 9999});
        requireNumberedEntry(cache, i);
    }
    assert(!cache.usesFlatStorage());
    {
        const auto last = Cache::flatPromotionThreshold - 1;
        const EdgeMask mask = numberedMask(last);
        cache.emplace(mask, {static_cast<int>(last), 7});
    }
    assert(cache.usesFlatStorage());
    assert(cache.size() == Cache::flatPromotionThreshold);
    assert(cache.flatCapacity() == 1024);
    assert(cache.overflowSize() == 0);
    // Promotion releases the previous map's nodes and bucket array.
    assert(cache.retainedBytes() == cache.flatCapacity() * 24);
    for (std::uint64_t i = 0; i < Cache::flatPromotionThreshold; ++i)
        requireNumberedEntry(cache, i);
    for (std::uint64_t i = Cache::flatPromotionThreshold; i < 1000; ++i)
    {
        const EdgeMask mask = numberedMask(i);
        cache.emplace(mask, {static_cast<int>(i), 7});
    }
    const auto retainedCapacity = cache.flatCapacity();
    assert(retainedCapacity == 2048);
    cache.clear();
    assert(!cache.usesFlatStorage());
    assert(cache.size() == 0);
    assert(cache.flatCapacity() == retainedCapacity);
    assert(cache.find(numberedMask(511)) == cache.end());

    for (std::uint64_t i = 2000; i < 2000 + Cache::flatPromotionThreshold; ++i)
    {
        assert(!cache.usesFlatStorage());
        const EdgeMask mask = numberedMask(i);
        cache.emplace(mask, {static_cast<int>(i), 7});
    }
    assert(cache.usesFlatStorage());
    assert(cache.flatCapacity() == retainedCapacity);
    assert(cache.overflowSize() == 0);
    assert(cache.retainedBytes() == retainedCapacity * 24);
    assert(cache.find(numberedMask(511)) == cache.end());
    for (std::uint64_t i = 2000; i < 2000 + Cache::flatPromotionThreshold; ++i)
        requireNumberedEntry(cache, i);

    // Retained flat slots must never activate a later, wider mask domain.
    cache.clear();
    EdgeMask::configure(129);
    for (std::uint64_t i = 0; i < 1024; ++i)
    {
        const std::array words{i & 255, (i >> 8) & 1, i >> 9};
        const EdgeMask mask = EdgeMask::fromActiveWords(words.data());
        cache.emplace(mask, {static_cast<int>(i), 7});
    }
    assert(!cache.usesFlatStorage());
    assert(cache.size() == 1024);
    assert(cache.overflowSize() == 1024);
    assert(cache.flatCapacity() == retainedCapacity);
    for (std::uint64_t i = 0; i < 1024; ++i)
    {
        const std::array words{i & 255, (i >> 8) & 1, i >> 9};
        const EdgeMask mask = EdgeMask::fromActiveWords(words.data());
        const auto found = cache.find(mask);
        assert(found != cache.end());
        assert(found->second.first == static_cast<int>(i));
    }
    cache.clear();
}

void testGrowthAndOverflow()
{
    EdgeMask::configure(128);
    Cache cache(Cache::policy::flat);
    constexpr std::uint64_t count = 80000;
    for (std::uint64_t i = 0; i < count; ++i)
    {
        // Repeated low words explicitly exercise second-word identity.
        const std::array words{i & 255, i >> 8};
        EdgeMask mask = EdgeMask::fromActiveWords(words.data());
        cache.emplace(mask, {static_cast<int>(i), 1});
    }
    assert(cache.size() == count);
    assert(cache.usesFlatStorage());
    assert(cache.flatCapacity() == Cache::maximumFlatCapacity);
    assert(cache.overflowSize() == count - Cache::maximumFlatCapacity * 7 / 10);
    const auto retained = cache.retainedBytes();
    for (std::uint64_t i = count; i-- > 0;)
    {
        const std::array words{i & 255, i >> 8};
        EdgeMask mask = EdgeMask::fromActiveWords(words.data());
        const auto found = cache.find(mask);
        assert(found != cache.end());
        assert(found->second.first == static_cast<int>(i));
        cache.emplace(mask, {1234, 9});
        const auto repeated = cache.find(mask);
        assert(repeated != cache.end());
        assert(repeated->second.first == static_cast<int>(i));
    }
    assert(cache.size() == count);
    assert(cache.retainedBytes() == retained);
    cache.clear();
    assert(cache.empty());
    assert(!cache.usesFlatStorage());
    assert(cache.retainedBytes() < retained);
}

int main()
{
    for (const auto policy : {Cache::policy::unordered, Cache::policy::flat})
        for (std::size_t width : {0, 1, 63, 64, 65, 127, 128, 129, 257})
            testWidth(width, policy);
    testPromotionAndGenerationReset(64);
    testPromotionAndGenerationReset(128);
    testGrowthAndOverflow();
}
