#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

/** A root occurrence stored without an EdgeMask from the producer thread. */
struct rootOccurrenceDescriptor
{
    std::size_t wordOffset = 0;
    std::int32_t fragment = 0;
};

/** One independently searchable root matching, expressed only as indices. */
struct rootJobDescriptor
{
    std::size_t firstOccurrence = 0;
    std::size_t secondOccurrence = 0;
    std::int32_t canonicalId = unknownCanonicalId;
    std::uint32_t duplicateSize = 0;
};

static_assert(std::is_trivially_copyable_v<rootOccurrenceDescriptor>);
static_assert(std::is_trivially_copyable_v<rootJobDescriptor>);

/**
 * Mask-free canonical state captured after the one root enumeration.
 *
 * Workers copy this seed into their local canonical caches. The tree interner
 * must travel with graphHashes because tree graph keys contain interner IDs.
 */
struct canonicalisationSeed
{
    decltype(graphHashMap) graphHashes;
    decltype(treeCanonAtomInterner) atomInterner;
    decltype(treeCanonLeafInterner) leafInterner;
    decltype(treeCanonInterner) treeInterner;
};

/**
 * Immutable input shared by every worker after construction.
 *
 * In particular, this type deliberately contains no EdgeMask, assemblyState,
 * or assemblyFragment. Wide masks belong to a thread-local arena and must be
 * reconstructed from occurrenceWords by the worker that will destroy them.
 */
struct SearchContext
{
    molGraph processedMolecule;
    vector<edgeL> universeEdges;
    vector<dagLevel> dag;
    vector<rootOccurrenceDescriptor> rootOccurrences;
    vector<std::uint64_t> occurrenceWords;
    vector<rootJobDescriptor> rootJobs;
    vector<int> homogeneousPathEdgePositions;
    canonicalisationSeed canonicalSeed;
    clock_t startedAt = 0;
    unsigned int bondCount = 0;
    int componentCount = 1;
    int rootAssemblyIndex = -1;
    bool enumerationLimit = false;

    SearchContext() = default;
    SearchContext(const SearchContext &) = delete;
    SearchContext &operator=(const SearchContext &) = delete;
    SearchContext(SearchContext &&) = default;
    SearchContext &operator=(SearchContext &&) = default;

    [[nodiscard]] std::size_t maskWordCount() const noexcept
    {
        return universeEdges.size() / EdgeMask::wordBits +
            (universeEdges.size() % EdgeMask::wordBits != 0);
    }
};

struct assemblyPathWitness
{
    vector<assemblyPathStep> current;
    vector<assemblyPathStep> best;
};

/** Mutable caches and fragmentation scratch owned by exactly one worker. */
struct assemblySearchStorage
{
    assemblyTranspositionTable states;
    assemblyPathWitness *pathway = nullptr;
    duplicateClassIndexWorkspace duplicateClassIndex;
    vi candidateKey;

    explicit assemblySearchStorage(assemblyPathWitness *_pathway = nullptr):
        states(1024),
        pathway(_pathway)
    {
        // One search-wide scratch key is safe because every table operation
        // finishes before the synchronous recursive call can reuse it.
        candidateKey.reserve(univEdgeList.size() + 1);
        if (pathway != nullptr)
        {
            pathway->current.reserve(univEdgeList.size());
            pathway->best.reserve(univEdgeList.size());
        }
    }
};

/**
 * Per-worker mutable search state. Construct only after configuring that
 * worker's EdgeMask and AtomMask domains from the shared SearchContext.
 */
struct WorkerContext
{
    ufdsMaskWorkspace fragmentation;
    assemblySearchStorage search;
    assemblyState root;
    assemblyState candidate;
    int assemblyIndex;

    explicit WorkerContext(const SearchContext &context):
        fragmentation(
            context.processedMolecule.mg.size(),
            context.universeEdges.size()
        ),
        assemblyIndex(context.rootAssemblyIndex)
    {
        fragmentation.homogeneousPathEdgePositions =
            context.homogeneousPathEdgePositions;

        EdgeMask rootMask;
        rootMask.set();
        root.appendFragment(
            rootMask,
            static_cast<int>(context.universeEdges.size()),
            unknownCanonicalId,
            false
        );
        root.assemblyHashCalculator(search.candidateKey);
        static_cast<void>(search.states.consider(search.candidateKey, 0));
        candidate.reserveFragments(3);
    }

    WorkerContext(const WorkerContext &) = delete;
    WorkerContext &operator=(const WorkerContext &) = delete;
    WorkerContext(WorkerContext &&) = delete;
    WorkerContext &operator=(WorkerContext &&) = delete;
};
