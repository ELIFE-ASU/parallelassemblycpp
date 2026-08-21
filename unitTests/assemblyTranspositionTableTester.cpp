// This executable is also built by Release/CI presets; keep its checks active.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../v5/assemblyTranspositionTable.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory_resource>
#include <random>
#include <unordered_map>
#include <vector>

using tableResult = assemblyTranspositionTable::result;

struct vectorHash
{
    std::size_t operator()(const std::vector<int> &key) const noexcept
    {
        std::size_t result = key.size();
        for (const int value : key)
        {
            result ^= static_cast<std::size_t>(
                static_cast<std::uint32_t>(value)
            ) + 0x9e3779b9 + (result << 6) + (result >> 2);
        }
        return result;
    }
};

class countingMemoryResource final : public std::pmr::memory_resource
{
public:
    std::size_t allocationCalls = 0;
    std::size_t deallocationCalls = 0;
    std::size_t allocatedBytes = 0;

private:
    void *do_allocate(std::size_t bytes, std::size_t alignment) override
    {
        ++allocationCalls;
        allocatedBytes += bytes;
        return std::pmr::new_delete_resource()->allocate(bytes, alignment);
    }

    void do_deallocate(
        void *memory,
        std::size_t bytes,
        std::size_t alignment
    ) override
    {
        ++deallocationCalls;
        std::pmr::new_delete_resource()->deallocate(
            memory,
            bytes,
            alignment
        );
    }

    bool do_is_equal(
        const std::pmr::memory_resource &other
    ) const noexcept override
    {
        return this == &other;
    }
};

void testBasicResultsAndExactKeys()
{
    assemblyTranspositionTable table(8);
    assert(table.size() == 0);
    assert(table.capacity() == 8);

    const std::array<int, 3> key{4, 1, 2};
    assert(table.consider(key, 7) == tableResult::inserted);
    assert(table.consider(key, 7) == tableResult::dominated);
    assert(table.consider(key, 6) == tableResult::dominated);
    assert(table.consider(key, 8) == tableResult::improved);
    assert(table.consider(key, 8) == tableResult::dominated);

    // The distinguished first value and the complete key length participate
    // in equality; neither a different head nor a same-prefix key aliases.
    const std::array<int, 3> differentHead{5, 1, 2};
    const std::array<int, 2> prefix{4, 1};
    const std::array<int, 4> extension{4, 1, 2, 0};
    assert(table.consider(differentHead, 1) == tableResult::inserted);
    assert(table.consider(prefix, 1) == tableResult::inserted);
    assert(table.consider(extension, 1) == tableResult::inserted);

    const std::span<const int> empty;
    assert(table.consider(empty, 0) == tableResult::inserted);
    assert(table.consider(empty, -1) == tableResult::dominated);
    assert(table.size() == 5);
}

void testScratchCopyAndHitAllocations()
{
    countingMemoryResource resource;
    assemblyTranspositionTable table(8, &resource);

    std::vector<int> scratch(4096);
    for (std::size_t i = 0; i < scratch.size(); i++)
        scratch[i] = static_cast<int>(i * 17 + 3);
    const std::vector<int> original = scratch;

    assert(table.consider(scratch, 10) == tableResult::inserted);
    assert(resource.allocationCalls > 0);
    const std::size_t allocationsAfterMiss = resource.allocationCalls;
    const std::size_t bytesAfterMiss = resource.allocatedBytes;

    // If the table retained the borrowed span, changing scratch would change
    // its stored key and the original value below would become another miss.
    std::fill(scratch.begin(), scratch.end(), -99);
    assert(table.consider(original, 10) == tableResult::dominated);
    assert(table.consider(original, 11) == tableResult::improved);
    assert(table.consider(original, 11) == tableResult::dominated);
    assert(table.size() == 1);

    // Dominated and improved hits must not request additional arena storage.
    assert(resource.allocationCalls == allocationsAfterMiss);
    assert(resource.allocatedBytes == bytesAfterMiss);
}

void testGrowthPreservesEntries()
{
    assemblyTranspositionTable table(8);
    constexpr int entryCount = 12000;
    for (int index = 0; index < entryCount; index++)
    {
        const std::array<int, 4> key{
            index % 23,
            index,
            index * 3,
            -index
        };
        assert(table.consider(key, index % 101) == tableResult::inserted);
    }
    assert(table.size() == entryCount);
    assert(table.capacity() > 8);

    for (int index = 0; index < entryCount; index++)
    {
        const std::array<int, 4> key{
            index % 23,
            index,
            index * 3,
            -index
        };
        assert(
            table.consider(key, index % 101) == tableResult::dominated
        );
    }
    assert(table.size() == entryCount);
}

void testRandomisedDifferential()
{
    assemblyTranspositionTable table(8);
    std::unordered_map<std::vector<int>, int, vectorHash> reference;
    std::vector<std::vector<int>> knownKeys;
    std::mt19937_64 random(0x6f70656e41646472ULL);

    constexpr std::size_t operationCount = 60000;
    for (std::size_t operation = 0; operation < operationCount; operation++)
    {
        std::vector<int> scratch;
        if (!knownKeys.empty() && random() % 100 < 67)
        {
            scratch = knownKeys[random() % knownKeys.size()];
        }
        else
        {
            const std::size_t length = random() % 25;
            scratch.resize(length);
            for (int &value : scratch)
                value = static_cast<int>(random() % 401) - 200;
        }
        const int score = static_cast<int>(random() % 301) - 100;

        auto [entry, inserted] = reference.try_emplace(scratch, score);
        tableResult expected;
        if (inserted)
        {
            expected = tableResult::inserted;
            knownKeys.push_back(scratch);
        }
        else if (score > entry->second)
        {
            entry->second = score;
            expected = tableResult::improved;
        }
        else expected = tableResult::dominated;

        assert(table.consider(scratch, score) == expected);
        assert(table.size() == reference.size());

        // Every operation uses disposable scratch storage. In particular, a
        // newly inserted key must remain intact after this mutation.
        for (int &value : scratch) value ^= 0x55aa55aa;
    }

    assert(table.capacity() > 8);
    for (auto &[key, bestScore] : reference)
    {
        assert(table.consider(key, bestScore) == tableResult::dominated);
        assert(table.consider(key, bestScore + 1) == tableResult::improved);
        assert(
            table.consider(key, bestScore + 1) == tableResult::dominated
        );
    }
    assert(table.size() == reference.size());
}

int main()
{
    testBasicResultsAndExactKeys();
    testScratchCopyAndHitAllocations();
    testGrowthPreservesEntries();
    testRandomisedDifferential();
}
