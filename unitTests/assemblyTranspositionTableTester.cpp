// This executable is also built by Release/CI presets; keep its checks active.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../v5/assemblyTranspositionTable.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <barrier>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory_resource>
#include <random>
#include <thread>
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

void testSharedExactKeysAndHashCollisions()
{
    sharedAssemblyTranspositionTable table;

    // These distinct keys deliberately have the same complete 32-bit hash,
    // and therefore enter the same shard and probe sequence. Equality must
    // still compare the borrowed key's full contents and length.
    const std::array<int, 4> firstCollision{
        -1744324134,
        -1879786136,
        873751343,
        1729211343
    };
    const std::array<int, 3> secondCollision{
        1933699411,
        -1276699930,
        -106575768
    };
    assert(
        assemblyTranspositionTable::keyHash(firstCollision) ==
        assemblyTranspositionTable::keyHash(secondCollision)
    );
    assert(table.consider(firstCollision, 4) == tableResult::inserted);
    assert(table.consider(secondCollision, 9) == tableResult::inserted);
    assert(table.consider(firstCollision, 4) == tableResult::dominated);
    assert(table.consider(secondCollision, 9) == tableResult::dominated);
    assert(table.consider(firstCollision, 5) == tableResult::improved);

    const std::array<int, 3> key{4, 1, 2};
    const std::array<int, 3> differentHead{5, 1, 2};
    const std::array<int, 2> prefix{4, 1};
    const std::array<int, 4> extension{4, 1, 2, 0};
    assert(table.consider(key, 1) == tableResult::inserted);
    assert(table.consider(differentHead, 1) == tableResult::inserted);
    assert(table.consider(prefix, 1) == tableResult::inserted);
    assert(table.consider(extension, 1) == tableResult::inserted);
    assert(table.size() == 6);
}

void testSharedBorrowedScratchLifetime()
{
    sharedAssemblyTranspositionTable table;
    std::vector<int> scratch(4096);
    for (std::size_t index = 0; index < scratch.size(); ++index)
        scratch[index] = static_cast<int>(index * 29 + 11);
    const std::vector<int> original = scratch;

    assert(table.consider(scratch, 12) == tableResult::inserted);
    std::fill(scratch.begin(), scratch.end(), -77);

    // A shared-table miss must copy the key before the shard lock is released;
    // it cannot retain worker scratch that will be reused immediately.
    assert(table.consider(original, 12) == tableResult::dominated);
    assert(table.consider(original, 13) == tableResult::improved);
    assert(table.consider(original, 13) == tableResult::dominated);
    assert(table.size() == 1);
}

void testSharedBestScoreSupportsLocalPromotion()
{
    sharedAssemblyTranspositionTable shared;
    assemblyTranspositionTable workerLocal(8);
    const std::array<int, 4> key{19, 3, 5, 7};

    auto sharedResult = shared.considerWithBest(key, 10);
    assert(sharedResult.outcome == tableResult::inserted);
    assert(sharedResult.bestSumDupBonds == 10);

    // A second worker reaches the same state with a weaker score. Its L1 miss
    // must be checked in L2, which returns the exact stronger score so the
    // worker can promote its local entry before pruning.
    assert(workerLocal.consider(key, 7) == tableResult::inserted);
    sharedResult = shared.considerWithBest(key, 7);
    assert(sharedResult.outcome == tableResult::dominated);
    assert(sharedResult.bestSumDupBonds == 10);
    assert(
        workerLocal.consider(key, sharedResult.bestSumDupBonds) ==
        tableResult::improved
    );
    assert(workerLocal.consider(key, 9) == tableResult::dominated);

    sharedResult = shared.considerWithBest(key, 12);
    assert(sharedResult.outcome == tableResult::improved);
    assert(sharedResult.bestSumDupBonds == 12);
    sharedResult = shared.considerWithBest(key, 11);
    assert(sharedResult.outcome == tableResult::dominated);
    assert(sharedResult.bestSumDupBonds == 12);
    assert(shared.size() == 1);
}

void testSharedConcurrentSameKeyMonotonicUpdate()
{
    constexpr int threadCount = 8;
    constexpr int rounds = 2048;
    sharedAssemblyTranspositionTable table(threadCount);
    const std::array<int, 6> key{7, 2, 3, 5, 11, 13};
    std::barrier start(threadCount);
    std::atomic<int> insertionCount{0};
    std::vector<std::thread> workers;
    workers.reserve(threadCount);

    for (int worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.arrive_and_wait();
            for (int round = 0; round < rounds; ++round)
            {
                const int score = round * threadCount + worker;
                if (
                    table.considerWithBestForWorker(
                        key,
                        score,
                        static_cast<std::size_t>(worker)
                    ).outcome == tableResult::inserted
                )
                    insertionCount.fetch_add(1, std::memory_order_relaxed);
            }
        });
    }
    for (std::thread &worker : workers) worker.join();

    const int maximumScore = rounds * threadCount - 1;
    assert(insertionCount.load(std::memory_order_relaxed) == 1);
    assert(table.size() == 1);
    assert(
        table.considerWithBestForWorker(key, maximumScore, 0).outcome ==
        tableResult::dominated
    );
    assert(
        table.considerWithBestForWorker(key, maximumScore + 1, 0).outcome ==
        tableResult::improved
    );
    assert(
        table.considerWithBestForWorker(key, maximumScore + 1, 0).outcome ==
        tableResult::dominated
    );
}

void testSharedIndependentKeyShardStress()
{
    sharedAssemblyTranspositionTable table;
    constexpr int threadCount = 8;
    constexpr int keysPerThread = 1024;
    std::array<bool, sharedAssemblyTranspositionTable::shardCount>
        visitedShards{};

    auto keyFor = [](int worker, int index)
    {
        const int ordinal = worker * keysPerThread + index;
        return std::array<int, 5>{
            0x13579bdf,
            worker,
            index,
            ordinal,
            -ordinal - 1
        };
    };
    for (int worker = 0; worker < threadCount; ++worker)
    {
        for (int index = 0; index < keysPerThread; ++index)
        {
            const auto key = keyFor(worker, index);
            const std::size_t shard =
                assemblyTranspositionTable::keyHash(key) &
                (sharedAssemblyTranspositionTable::shardCount - 1);
            visitedShards[shard] = true;
        }
    }
    assert(std::all_of(
        visitedShards.begin(),
        visitedShards.end(),
        [](bool visited) {return visited;}
    ));

    std::barrier start(threadCount);
    std::atomic<int> failureCount{0};
    std::vector<std::thread> workers;
    workers.reserve(threadCount);
    for (int worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.arrive_and_wait();
            for (int index = 0; index < keysPerThread; ++index)
            {
                const auto key = keyFor(worker, index);
                const int score = worker * keysPerThread + index;
                if (table.consider(key, score) != tableResult::inserted)
                    failureCount.fetch_add(1, std::memory_order_relaxed);
                if (table.consider(key, score) != tableResult::dominated)
                    failureCount.fetch_add(1, std::memory_order_relaxed);
                if (table.consider(key, score + 1) != tableResult::improved)
                    failureCount.fetch_add(1, std::memory_order_relaxed);
                if (table.consider(key, score + 1) != tableResult::dominated)
                    failureCount.fetch_add(1, std::memory_order_relaxed);
            }
        });
    }
    for (std::thread &worker : workers) worker.join();

    assert(failureCount.load(std::memory_order_relaxed) == 0);
    assert(table.size() == threadCount * keysPerThread);
    for (int worker = 0; worker < threadCount; ++worker)
    {
        for (int index = 0; index < keysPerThread; ++index)
        {
            const auto key = keyFor(worker, index);
            const int score = worker * keysPerThread + index + 1;
            assert(table.consider(key, score) == tableResult::dominated);
        }
    }
}

int main()
{
    testBasicResultsAndExactKeys();
    testScratchCopyAndHitAllocations();
    testGrowthPreservesEntries();
    testRandomisedDifferential();
    testSharedExactKeysAndHashCollisions();
    testSharedBorrowedScratchLifetime();
    testSharedBestScoreSupportsLocalPromotion();
    testSharedConcurrentSameKeyMonotonicUpdate();
    testSharedIndependentKeyShardStress();
}
