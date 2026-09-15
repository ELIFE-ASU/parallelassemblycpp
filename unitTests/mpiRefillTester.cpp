#define PARALLELASSEMBLYCPP_NO_MAIN
#include "../src/main.cpp"

#include <cassert>

namespace
{
constexpr int requestPostedTag = 19101;
constexpr int replyObservedTag = 19102;

ParallelSearchWorkerTelemetry snapshot(
    const MpiDistributedSearchController &controller
)
{
    ParallelSearchWorkerTelemetry telemetry;
    controller.captureTelemetry(telemetry);
    return telemetry;
}

void sendSignal(int destination, int tag)
{
    MPI_Send(nullptr, 0, MPI_BYTE, destination, tag, MPI_COMM_WORLD);
}

void receiveSignal(int source, int tag)
{
    MPI_Recv(
        nullptr, 0, MPI_BYTE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
}

void recordLease(vector<int> &visits, size_t begin, size_t end)
{
    assert(begin < end && end <= visits.size());
    for (size_t index = begin; index < end; ++index) ++visits[index];
}

void requireExactCoverage(const vector<int> &localVisits)
{
    vector<int> globalVisits(localVisits.size());
    MPI_Allreduce(
        localVisits.data(), globalVisits.data(),
        static_cast<int>(localVisits.size()), MPI_INT, MPI_SUM, MPI_COMM_WORLD
    );
    for (int visits : globalVisits) assert(visits == 1);
}

void requireDrained(const MpiDistributedSearchController &controller)
{
    const auto telemetry = snapshot(controller);
    assert(telemetry.mpiRefillRequests == telemetry.mpiRefillReplies);
    assert(telemetry.mpiPendingRefillsHighWatermark <= 1);
    assert(telemetry.mpiPrefetchedRefills <= telemetry.mpiRefillRequests);
    assert(telemetry.mpiProgressCalls > 0);
    assert(telemetry.mpiMaximumProgressGapNanoseconds > 0);
}

void testLowWatermarkAndReplyIncumbent()
{
    constexpr size_t totalJobs = 128;
    std::atomic<int> best{100};
    MpiDistributedSearchController controller(
        totalJobs, 8, 1, parallelAssemblyCppMpiRank, 2, best
    );
    vector<int> visits(totalJobs);
    controller.progress();
    std::this_thread::sleep_for(std::chrono::milliseconds(3));
    controller.progress();
    MPI_Barrier(MPI_COMM_WORLD);
    assert(snapshot(controller).mpiRefillRequests == 0);

    if (parallelAssemblyCppMpiRank == 0)
    {
        // MPI itself can progress while rank zero deliberately leaves the
        // refill broker unserviced. A blocking refill cannot reach this send.
        receiveSignal(1, requestPostedTag);
        best.store(7, std::memory_order_relaxed);
        int observed = 0;
        while (observed == 0)
        {
            controller.progress();
            MPI_Iprobe(
                1, replyObservedTag, MPI_COMM_WORLD, &observed, MPI_STATUS_IGNORE
            );
        }
        receiveSignal(1, replyObservedTag);
    }
    else
    {
        size_t begin = 0;
        size_t end = 0;
        assert(controller.claimRootLease(4, begin, end) ==
            distributedRootAvailability::lease);
        assert(begin == 8 && end == 12);
        recordLease(visits, begin, end);
        for (int iteration = 0; iteration < 4; ++iteration) controller.progress();
        const auto pending = snapshot(controller);
        assert(pending.mpiRefillRequests == 1);
        assert(pending.mpiPrefetchedRefills == 1);
        assert(pending.mpiRefillReplies == 0);
        assert(pending.mpiPendingRefillsHighWatermark == 1);
        assert(best.load(std::memory_order_relaxed) == 100);

        // Existing work remains claimable while its one refill is pending.
        assert(controller.claimRootLease(2, begin, end) ==
            distributedRootAvailability::lease);
        assert(begin == 12 && end == 14);
        recordLease(visits, begin, end);
        sendSignal(0, requestPostedTag);
        while (best.load(std::memory_order_relaxed) != 7) controller.progress();
        // The incumbent is usable immediately, but publishing the new range
        // before the old one is consumed would lose the final two old roots.
        assert(snapshot(controller).mpiRefillReplies == 0);
        assert(controller.claimRootLease(2, begin, end) ==
            distributedRootAvailability::lease);
        assert(begin == 14 && end == 16);
        recordLease(visits, begin, end);
        // Reaching the old end while a reply is staged is still not completion.
        assert(controller.claimRootLease(1, begin, end) ==
            distributedRootAvailability::wait);
        while (snapshot(controller).mpiRefillReplies == 0) controller.progress();
        assert(best.load(std::memory_order_relaxed) == 7);
        sendSignal(0, replyObservedTag);
    }

    while (true)
    {
        size_t begin = 0;
        size_t end = 0;
        const auto availability = controller.claimRootLease(1, begin, end);
        if (availability == distributedRootAvailability::complete) break;
        if (availability == distributedRootAvailability::lease)
            recordLease(visits, begin, end);
        controller.progress();
    }
    controller.waitForGlobalCompletion();
    requireDrained(controller);
    assert(snapshot(controller).mpiMaximumProgressGapNanoseconds >= 2000000);
    requireExactCoverage(visits);
}

void testCancellationWithPendingReply()
{
    std::atomic<int> best{100};
    MpiDistributedSearchController controller(
        128, 8, 1, parallelAssemblyCppMpiRank, 2, best
    );
    controller.progress();
    MPI_Barrier(MPI_COMM_WORLD);
    if (parallelAssemblyCppMpiRank == 0)
    {
        receiveSignal(1, requestPostedTag);
    }
    else
    {
        size_t begin = 0;
        size_t end = 0;
        assert(controller.claimRootLease(4, begin, end) ==
            distributedRootAvailability::lease);
        controller.progress();
        const auto pending = snapshot(controller);
        assert(pending.mpiRefillRequests == 1);
        assert(pending.mpiRefillReplies == 0);
        searchCancellationFlag.store(true, std::memory_order_release);
        sendSignal(0, requestPostedTag);
    }
    // Both ranks finish local work while rank one's receive is outstanding.
    // Rank zero must keep servicing it and keep reply buffers alive until the
    // response drains, including when the cancellation makes that reply empty.
    controller.waitForGlobalCompletion();
    requireDrained(controller);
    if (parallelAssemblyCppMpiRank == 1)
        assert(searchCancellationFlag.load(std::memory_order_acquire));
}

void testConcurrentClaimsAndPaddedTail()
{
    constexpr int localWorkers = 4;
    constexpr size_t leaseSize = 8;
    constexpr size_t totalJobs = 259;
    constexpr size_t initialBlock = leaseSize * localWorkers * 2;
    constexpr size_t paddedJobs =
        (totalJobs + initialBlock - 1) / initialBlock * initialBlock;
    std::atomic<int> best{100};
    MpiDistributedSearchController controller(
        totalJobs, leaseSize, localWorkers,
        parallelAssemblyCppMpiRank * localWorkers, localWorkers * 2, best
    );
    std::array<vector<int>, localWorkers> visits;
    for (auto &workerVisits : visits) workerVisits.resize(paddedJobs);
    std::atomic<int> finishedWorkers{0};
    std::array<std::thread, localWorkers - 1> workers;
    for (size_t index = 0; index < workers.size(); ++index)
    {
        workers[index] = std::thread([&, index]
        {
            while (true)
            {
                size_t begin = 0;
                size_t end = 0;
                const auto availability = controller.claimRootLease(
                    index + 1, begin, end
                );
                if (availability == distributedRootAvailability::complete) break;
                if (availability == distributedRootAvailability::lease)
                    recordLease(visits[index + 1], begin, end);
                else std::this_thread::yield();
            }
            finishedWorkers.fetch_add(1, std::memory_order_release);
        });
    }
    bool mainWorkerComplete = false;
    while (
        !mainWorkerComplete ||
        finishedWorkers.load(std::memory_order_acquire) != localWorkers - 1
    )
    {
        controller.progress();
        size_t begin = 0;
        size_t end = 0;
        const auto availability = controller.claimRootLease(4, begin, end);
        mainWorkerComplete = availability == distributedRootAvailability::complete;
        if (availability == distributedRootAvailability::lease)
            recordLease(visits[0], begin, end);
    }
    for (auto &worker : workers) worker.join();
    controller.waitForGlobalCompletion();
    requireDrained(controller);
    for (size_t index = 1; index < visits.size(); ++index)
        for (size_t root = 0; root < paddedJobs; ++root)
            visits[0][root] += visits[index][root];
    requireExactCoverage(visits[0]);
}

void testCancellationWithStagedReply()
{
    std::atomic<int> best{100};
    MpiDistributedSearchController controller(
        128, 8, 1, parallelAssemblyCppMpiRank, 2, best
    );
    controller.progress();
    MPI_Barrier(MPI_COMM_WORLD);
    if (parallelAssemblyCppMpiRank == 0)
    {
        receiveSignal(1, requestPostedTag);
        best.store(7, std::memory_order_relaxed);
        int observed = 0;
        while (observed == 0)
        {
            controller.progress();
            MPI_Iprobe(
                1, replyObservedTag, MPI_COMM_WORLD, &observed, MPI_STATUS_IGNORE
            );
        }
        receiveSignal(1, replyObservedTag);
    }
    else
    {
        size_t begin = 0;
        size_t end = 0;
        assert(controller.claimRootLease(4, begin, end) ==
            distributedRootAvailability::lease);
        controller.progress();
        sendSignal(0, requestPostedTag);
        while (best.load(std::memory_order_relaxed) != 7) controller.progress();
        // A nonempty reply has arrived, but four original slots are still
        // active. Cancellation must retire this staged reply without waiting
        // for workers to consume those original slots or publishing new work.
        assert(snapshot(controller).mpiRefillRequests == 1);
        assert(snapshot(controller).mpiRefillReplies == 0);
        searchCancellationFlag.store(true, std::memory_order_release);
        sendSignal(0, replyObservedTag);
    }
    controller.waitForGlobalCompletion();
    requireDrained(controller);
    if (parallelAssemblyCppMpiRank == 1)
    {
        assert(snapshot(controller).mpiRefillRequests == 1);
        size_t begin = 0;
        size_t end = 0;
        assert(controller.claimRootLease(4, begin, end) ==
            distributedRootAvailability::lease);
        assert(begin == 12 && end == 16);
    }
}

void testEmptyFrontier()
{
    std::atomic<int> best{0};
    MpiDistributedSearchController controller(
        0, 8, 1, parallelAssemblyCppMpiRank, 2, best
    );
    size_t begin = 0;
    size_t end = 0;
    assert(controller.claimRootLease(1, begin, end) ==
        distributedRootAvailability::complete);
    controller.waitForGlobalCompletion();
    requireDrained(controller);
    assert(snapshot(controller).mpiRefillRequests == 0);
}
}

int main(int argc, char **argv)
{
    int provided = MPI_THREAD_SINGLE;
    const int initialized = MPI_Init_thread(
        &argc, &argv, MPI_THREAD_FUNNELED, &provided
    );
    assert(initialized == MPI_SUCCESS && provided >= MPI_THREAD_FUNNELED);
    MPI_Comm_rank(MPI_COMM_WORLD, &parallelAssemblyCppMpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &parallelAssemblyCppMpiSize);
    assert(parallelAssemblyCppMpiSize == 2);
    searchTelemetryEnabled = true;
    testLowWatermarkAndReplyIncumbent();
    testCancellationWithPendingReply();
    searchCancellationFlag.store(false, std::memory_order_release);
    MPI_Barrier(MPI_COMM_WORLD);
    testConcurrentClaimsAndPaddedTail();
    testCancellationWithStagedReply();
    searchCancellationFlag.store(false, std::memory_order_release);
    MPI_Barrier(MPI_COMM_WORLD);
    testEmptyFrontier();
    if (parallelAssemblyCppMpiRank == 0)
        cout << "PASS MPI asynchronous refill, incumbent replies, progress gaps, "
                "concurrent claims, and cancellation drain\n";
    MPI_Finalize();
}
