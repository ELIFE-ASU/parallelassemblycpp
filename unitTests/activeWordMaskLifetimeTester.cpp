// This executable is also built by Release/CI presets; keep its checks active.
#ifdef NDEBUG
#undef NDEBUG
#endif

#define PARALLELASSEMBLYCPP_SEARCH_LOCAL thread_local
#include "../src/activeWordMask.h"

#include <array>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <thread>
#include <vector>

namespace
{
// Only one worker runs at a time, and the main thread inspects these counters
// after joining it. Neither the tracker nor the mask workload allocates storage
// for bookkeeping, so the tracked allocations belong to the mask arenas.
std::array<void *, 64> trackedAllocations{};
std::size_t liveAllocationCount = 0;
std::size_t allocationCount = 0;
thread_local bool trackNewAllocations = false;

void recordAllocation(void *memory)
{
    if (!trackNewAllocations) return;
    for (void *&allocation : trackedAllocations)
    {
        if (allocation != nullptr) continue;
        allocation = memory;
        ++liveAllocationCount;
        ++allocationCount;
        return;
    }
    std::abort();
}

void recordDeallocation(void *memory) noexcept
{
    if (memory == nullptr) return;
    for (void *&allocation : trackedAllocations)
    {
        if (allocation != memory) continue;
        allocation = nullptr;
        --liveAllocationCount;
        return;
    }
}

EdgeMask &persistentMask()
{
    thread_local EdgeMask mask;
    return mask;
}

std::vector<AtomMask> &persistentMasks()
{
    thread_local std::vector<AtomMask> masks;
    return masks;
}

template<typename Mask>
void fillMultipleArenaBlocks()
{
    const std::size_t initialAllocationCount = allocationCount;
    std::array<Mask, 4096> masks;
    for (Mask &mask : masks)
    {
        mask.set(Mask::size() - 1);
        assert(mask.count() == 1);
    }
    assert(allocationCount >= initialAllocationCount + 2);
}

void workerCalculation()
{
    // Register this TLS destructor before configuring the arena. Reconstruct
    // the mask around each width change, as required by configure(), while
    // retaining the original TLS destructor registration order.
    EdgeMask &persistent = persistentMask();
    std::vector<AtomMask> &persistentAtoms = persistentMasks();
    persistentAtoms.reserve(2);
    std::destroy_at(&persistent);

    // Thread runtime allocations have already happened; only arena work is
    // tracked from here until the end of the calculation.
    trackNewAllocations = true;
    EdgeMask::configure(129);
    AtomMask::configure(257);
    std::construct_at(&persistent);
    fillMultipleArenaBlocks<EdgeMask>();
    fillMultipleArenaBlocks<AtomMask>();
    persistent.set(128);
    persistentAtoms.emplace_back().set(256);

    std::destroy_at(&persistent);
    persistentAtoms.clear();
    EdgeMask::configure(1);
    AtomMask::configure(1);
    assert(liveAllocationCount == 0);

    EdgeMask::configure(513);
    AtomMask::configure(1025);
    std::construct_at(&persistent);
    fillMultipleArenaBlocks<EdgeMask>();
    fillMultipleArenaBlocks<AtomMask>();
    persistent.set(512);
    persistentAtoms.emplace_back().set(1024);
    persistentAtoms.push_back(persistentAtoms.front());
    {
        EdgeMask copy = persistent;
        copy.set(0);
        assert(copy.count() == 2);
        assert(persistent.count() == 1);
    }
    assert(liveAllocationCount > 0);
    trackNewAllocations = false;
    // Do not reset either arena: worker exit must reclaim both domains, and
    // persistent mask destructors must be safe even if the arena was already
    // torn down.
}
}

void *operator new(std::size_t size)
{
    void *memory = std::malloc(size == 0 ? 1 : size);
    if (memory == nullptr) throw std::bad_alloc();
    recordAllocation(memory);
    return memory;
}

void operator delete(void *memory) noexcept
{
    recordDeallocation(memory);
    std::free(memory);
}

void operator delete(void *memory, std::size_t) noexcept
{
    ::operator delete(memory);
}

int main()
{
    for (int calculation = 0; calculation < 4; ++calculation)
    {
        std::thread worker(workerCalculation);
        worker.join();
        if (liveAllocationCount != 0)
        {
            std::fprintf(
                stderr,
                "Worker %d leaked %zu mask arena allocations\n",
                calculation,
                liveAllocationCount
            );
            return 1;
        }
    }
}
