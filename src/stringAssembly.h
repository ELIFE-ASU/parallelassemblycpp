#pragma once

/*
 * Adapted for ParallelAssemblyCpp from the string assembly implementation in
 * croningroup/public/assemblycpp-public commit 2a87948, authored by Stuart
 * Marshall from work by Ian Seet and Leroy Cronin. See README.md and
 * License.md for source details and licensing.
 */

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <ctime>
#include <exception>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(PARALLELASSEMBLYCPP_USE_OPENMP)
#include <omp.h>
#endif

#include "additionChainBounds.h"
#include "clockTicks.h"
#include "stringEncoding.h"
#include "stringRepair.h"

namespace parallelassemblycpp::detail::stringAssembly
{

/** A half-open interval of Unicode scalar values in the original string. */
struct Interval
{
    int offset = 0;
    int length = 0;

    bool operator==(const Interval &) const = default;
};

struct IntervalHash
{
    size_t operator()(const Interval &interval) const noexcept
    {
        size_t result = std::hash<int>{}(interval.offset);
        result ^= std::hash<int>{}(interval.length) +
            static_cast<size_t>(0x9e3779b97f4a7c15ULL) +
            (result << 6) + (result >> 2);
        return result;
    }
};

struct IntegerVectorHash
{
    size_t operator()(const std::vector<int> &values) const noexcept
    {
        size_t result = values.size();
        for (const int value : values)
        {
            result ^= std::hash<int>{}(value) +
                static_cast<size_t>(0x9e3779b97f4a7c15ULL) +
                (result << 6) + (result >> 2);
        }
        return result;
    }
};

struct PathwayStep
{
    Interval match;
    Interval duplicate;
};

using CancellationCheck = bool (*)();

struct Options
{
    bool acceptReversed = false;
    unsigned long long runtimeTicks =
        std::numeric_limits<unsigned long long>::max();
    CancellationCheck cancellationRequested = nullptr;
    // Local OpenMP workers; independent processes can partition root jobs
    // through shardIndex/shardCount using the same input and options.
    int threadCount = 1;
    size_t shardIndex = 0;
    size_t shardCount = 1;
    // A nonnegative target stops serial reconstruction at its first witness.
    int targetAssemblyIndex = -1;
    // Replay a completed parallel/sharded search serially for deterministic
    // pathway ordering. Witnesses are retained even when replay is disabled.
    bool reconstructPathway = true;
};

struct Result
{
    // The empty string has index -1; a single symbol has index zero.
    int assemblyIndex = -1;
    unsigned long long clockTicks = 0;
    bool runtimeLimitReached = false;
    bool interrupted = false;
    std::vector<PathwayStep> pathway;
};

namespace implementation
{

using stringEncoding::decodeInput;
using stringEncoding::writeJsonString;

/**
 * Turn the Re-Pair grammar into the exact solver's duplicate-removal witness.
 * Visit parents before children, retaining one occurrence of each product and
 * expanding only that occurrence. Expanding discarded copies would reuse text
 * that is no longer present. Unreachable productions need no construction.
 */
inline std::vector<PathwayStep> repairPathway(
    const stringRepair::Result &grammar,
    const std::function<bool()> &shouldStop = {}
)
{
    std::size_t work = 0;
    const auto poll = [&] {
        if (shouldStop && (work++ % 1024 == 0) && shouldStop())
            throw stringRepair::implementation::Cancelled{};
    };
    const std::size_t symbolCount = grammar.terminals.size() + grammar.rules.size();
    std::vector<int> lengths(symbolCount, 1);
    std::vector<std::vector<stringRepair::Fragment>> occurrences(symbolCount);
    for (const auto &rule : grammar.rules)
    {
        poll();
        lengths[static_cast<std::size_t>(rule.id)] = rule.length;
    }
    const auto append = [&](const stringRepair::Fragment &fragment) {
        if (fragment.length > 1)
            occurrences[static_cast<std::size_t>(fragment.symbol)].push_back(fragment);
    };
    for (const auto &fragment : grammar.residual)
    {
        poll();
        append(fragment);
    }

    std::vector<PathwayStep> pathway;
    for (auto rule = grammar.rules.rbegin(); rule != grammar.rules.rend(); ++rule)
    {
        poll();
        const auto &copies = occurrences[static_cast<std::size_t>(rule->id)];
        if (copies.empty()) continue;
        std::size_t representative = 0;
        for (std::size_t i = 1; i < copies.size(); ++i)
        {
            poll();
            if (copies[i].offset < copies[representative].offset) representative = i;
        }
        const auto &kept = copies[representative];
        for (std::size_t i = 0; i < copies.size(); ++i)
        {
            poll();
            if (i != representative)
                pathway.push_back({{kept.offset, kept.length},
                    {copies[i].offset, copies[i].length}});
        }

        // Reversing a production swaps its children as well as their senses.
        const int left = kept.reversed ? rule->right : rule->left;
        const int right = kept.reversed ? rule->left : rule->right;
        const bool leftReversed = kept.reversed !=
            (kept.reversed ? rule->rightReversed : rule->leftReversed);
        const bool rightReversed = kept.reversed !=
            (kept.reversed ? rule->leftReversed : rule->rightReversed);
        const int leftLength = lengths[static_cast<std::size_t>(left)];
        append({left, leftReversed, kept.offset, leftLength});
        append({right, rightReversed, kept.offset + leftLength,
            lengths[static_cast<std::size_t>(right)]});
    }
    return pathway;
}

/** Sorted union of intervals having one fixed insertion length. */
class FixedIntervalMap
{
    std::map<int, int> intervals_;
    int covered_ = 0;
    int insertionLength_ = 2;

public:
    FixedIntervalMap() = default;

    explicit FixedIntervalMap(int insertionLength):
        insertionLength_(insertionLength) {}

    void insert(int offset)
    {
        int mergedBegin = offset;
        int mergedEnd = offset + insertionLength_;
        auto next = intervals_.lower_bound(offset);

        if (next != intervals_.begin())
        {
            auto previous = std::prev(next);
            const int previousEnd = previous->first + previous->second;
            if (previousEnd >= mergedBegin)
            {
                mergedBegin = previous->first;
                mergedEnd = std::max(mergedEnd, previousEnd);
                covered_ -= previous->second;
                next = intervals_.erase(previous);
            }
        }

        while (next != intervals_.end() && next->first <= mergedEnd)
        {
            mergedEnd = std::max(mergedEnd, next->first + next->second);
            covered_ -= next->second;
            next = intervals_.erase(next);
        }

        const int mergedLength = mergedEnd - mergedBegin;
        intervals_.emplace(mergedBegin, mergedLength);
        covered_ += mergedLength;
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return intervals_.empty();
    }

    [[nodiscard]] bool containsTwoDisjointInsertions() const noexcept
    {
        return static_cast<long long>(covered_) >=
            2LL * static_cast<long long>(insertionLength_);
    }

    [[nodiscard]] std::vector<Interval> intervals(
        bool omitSingletons = false
    ) const
    {
        std::vector<Interval> result;
        result.reserve(intervals_.size());
        for (const auto &[offset, length] : intervals_)
        {
            if (!omitSingletons || length > 1)
                result.push_back({offset, length});
        }
        return result;
    }
};

struct PotentialDuplicate
{
    Interval interval;
    size_t fragment = 0;

    void extend(
        std::vector<PotentialDuplicate> &output,
        const Interval &containingFragment
    ) const
    {
        if (
            interval.offset + interval.length <
            containingFragment.offset + containingFragment.length
        )
        {
            PotentialDuplicate extended = *this;
            extended.interval.length++;
            output.push_back(extended);
        }
    }

    [[nodiscard]] bool overlaps(const PotentialDuplicate &other) const noexcept
    {
        const int end = interval.offset + interval.length;
        const int otherEnd = other.interval.offset + other.interval.length;
        return interval.offset < otherEnd && other.interval.offset < end;
    }
};

struct ValidMatching
{
    Interval first;
    Interval second;
    size_t firstFragment = 0;
    size_t secondFragment = 0;
    int fragmentLength = 0;
};

class DuplicateSet
{
    size_t fragmentLength_ = 0;
    std::vector<PotentialDuplicate> occurrences_;
    int minimumOffset_ = 0;
    int maximumOffset_ = 0;
    bool multipleFragments_ = false;

public:
    explicit DuplicateSet(size_t fragmentLength):
        fragmentLength_(fragmentLength) {}

    void insert(const PotentialDuplicate &occurrence)
    {
        if (occurrences_.empty())
        {
            minimumOffset_ = occurrence.interval.offset;
            maximumOffset_ = occurrence.interval.offset;
        }
        else
        {
            minimumOffset_ = std::min(minimumOffset_, occurrence.interval.offset);
            maximumOffset_ = std::max(maximumOffset_, occurrence.interval.offset);
            multipleFragments_ = multipleFragments_ ||
                occurrence.fragment != occurrences_.front().fragment;
        }
        occurrences_.push_back(occurrence);
    }

    [[nodiscard]] bool isValid() const noexcept
    {
        return !occurrences_.empty() &&
            (multipleFragments_ ||
                maximumOffset_ - minimumOffset_ >= static_cast<int>(fragmentLength_));
    }

    /** Visits pairs in reverse insertion order without storing the pairs. */
    class MatchingCursor
    {
        const DuplicateSet *duplicates_ = nullptr;
        // One past the first occurrence; second_ descends over its partners.
        size_t first_ = 0;
        size_t second_ = 0;

    public:
        MatchingCursor() = default;

        explicit MatchingCursor(const DuplicateSet &duplicates):
            duplicates_(&duplicates),
            first_(duplicates.isValid() ? duplicates.occurrences_.size() - 1 : 0),
            second_(duplicates.occurrences_.size()) {}

        bool next(ValidMatching &matching)
        {
            while (first_ > 0)
            {
                while (second_ > first_)
                {
                    const PotentialDuplicate &left =
                        duplicates_->occurrences_[first_ - 1];
                    const PotentialDuplicate &right =
                        duplicates_->occurrences_[--second_];
                    if (
                        left.fragment == right.fragment &&
                        left.overlaps(right)
                    ) continue;
                    matching = {
                        left.interval, right.interval,
                        left.fragment, right.fragment,
                        static_cast<int>(duplicates_->fragmentLength_)
                    };
                    return true;
                }
                first_--;
                second_ = duplicates_->occurrences_.size();
            }
            return false;
        }
    };

    /** The set must remain alive and unchanged while its cursor is in use. */
    [[nodiscard]] MatchingCursor matchingCursor() const
    {
        return MatchingCursor(*this);
    }

    bool extendValidOccurrences(
        std::vector<PotentialDuplicate> &output,
        const std::vector<Interval> &fragments,
        std::vector<FixedIntervalMap> *survivingIntervals = nullptr
    ) const
    {
        bool extendedAny = false;
        for (const PotentialDuplicate &occurrence : occurrences_)
        {
            // Different fragments always supply a partner. In one fragment,
            // an equal-length occurrence has a disjoint partner exactly when
            // it is far enough from one of the two extreme offsets.
            const int offset = occurrence.interval.offset;
            if (
                !multipleFragments_ &&
                offset - minimumOffset_ < static_cast<int>(fragmentLength_) &&
                maximumOffset_ - offset < static_cast<int>(fragmentLength_)
            ) continue;

            if (survivingIntervals != nullptr)
            {
                (*survivingIntervals)[occurrence.fragment].insert(offset);
            }
            occurrence.extend(output, fragments[occurrence.fragment]);
            extendedAny = true;
        }
        return extendedAny;
    }
};

struct AssemblyState
{
    std::vector<Interval> intervals;
    int duplicatedSymbols = 0;
};

struct Enumeration
{
    std::vector<std::map<int, DuplicateSet>> duplicateSetsByLength;
    std::vector<std::vector<Interval>> remnantIntervals;
};

/** Shared coordination only; canonical IDs and transposition tables are local. */
struct ParallelControl
{
    explicit ParallelControl(int initialIndex): bestAssemblyIndex(initialIndex) {}

    std::atomic<int> bestAssemblyIndex;
    std::atomic<bool> runtimeLimitReached{false};
    std::atomic<bool> interrupted{false};
    std::atomic<bool> failed{false};
    std::mutex cancellationMutex;
    std::mutex failureMutex;
    std::exception_ptr failure;
};

class Search
{
    std::u32string original_;
    Options options_;
    std::clock_t started_ = 0;
    std::unordered_map<std::u32string, int> stringIds_;
    std::unordered_map<Interval, int, IntervalHash> intervalIds_;
    std::unordered_map<int, std::vector<int>> rollingHash_;
    std::unordered_map<std::vector<int>, int, IntegerVectorHash> states_;
    std::vector<PathwayStep> currentPath_;
    std::vector<PathwayStep> bestPath_;
    int bestAssemblyIndex_ = -1;
    // A whole-input composition bound, never a sum over fragments: fragment
    // constructions can share products. Copies/replay reuse this same proof.
    int assemblyLowerBound_ = -1;
    bool runtimeLimitReached_ = false;
    bool interrupted_ = false;
    bool targetReached_ = false;
    ParallelControl *parallelControl_ = nullptr;

    [[nodiscard]] bool optimumReached() const noexcept
    {
        return assemblyLowerBound_ >= 0 &&
            (parallelControl_ == nullptr
                ? bestAssemblyIndex_
                : parallelControl_->bestAssemblyIndex.load(std::memory_order_relaxed))
                <= assemblyLowerBound_;
    }

    [[nodiscard]] unsigned long long elapsedTicks() const noexcept
    {
        return assembly_clock::budgetTicks(started_, std::clock());
    }

    bool shouldStop()
    {
        if (targetReached_) return true;
        if (optimumReached()) return true;
        if (parallelControl_ != nullptr)
        {
            ParallelControl &control = *parallelControl_;
            if (
                control.failed.load(std::memory_order_relaxed) ||
                control.interrupted.load(std::memory_order_relaxed) ||
                control.runtimeLimitReached.load(std::memory_order_relaxed)
            ) return true;
            if (
                options_.cancellationRequested == nullptr &&
                options_.runtimeTicks ==
                    std::numeric_limits<unsigned long long>::max()
            ) return false;

            // A cancellation callback need not itself be thread-safe. Only
            // one worker polls it at a time, and cancellation takes priority
            // over a runtime limit observed by the same poll.
            std::lock_guard<std::mutex> lock(control.cancellationMutex);
            if (
                control.failed.load(std::memory_order_relaxed) ||
                control.interrupted.load(std::memory_order_relaxed) ||
                control.runtimeLimitReached.load(std::memory_order_relaxed)
            ) return true;
            if (
                options_.cancellationRequested != nullptr &&
                options_.cancellationRequested()
            )
            {
                control.interrupted.store(true, std::memory_order_relaxed);
                return true;
            }
            if (
                options_.runtimeTicks !=
                    std::numeric_limits<unsigned long long>::max() &&
                elapsedTicks() >= options_.runtimeTicks
            )
            {
                control.runtimeLimitReached.store(true, std::memory_order_relaxed);
                return true;
            }
            return false;
        }
        if (interrupted_ || runtimeLimitReached_) return true;
        if (
            options_.cancellationRequested != nullptr &&
            options_.cancellationRequested()
        )
        {
            interrupted_ = true;
            return true;
        }
        if (
            options_.runtimeTicks !=
                std::numeric_limits<unsigned long long>::max() &&
            elapsedTicks() >= options_.runtimeTicks
        )
        {
            runtimeLimitReached_ = true;
            return true;
        }
        return false;
    }

    [[nodiscard]] int pruningIndex() const noexcept
    {
        int bound = parallelControl_ == nullptr
            ? bestAssemblyIndex_
            : parallelControl_->bestAssemblyIndex.load(std::memory_order_relaxed);
        if (options_.targetAssemblyIndex >= 0 && options_.targetAssemblyIndex < bound)
            bound = options_.targetAssemblyIndex + 1;
        return bound;
    }

    void recordBest(int assemblyIndex, const std::vector<PathwayStep> &pathway)
    {
        if (assemblyIndex < bestAssemblyIndex_)
        {
            // Keep the witness before publishing its bound to other workers.
            bestPath_ = pathway;
            bestAssemblyIndex_ = assemblyIndex;
            if (parallelControl_ != nullptr)
            {
                int previous = parallelControl_->bestAssemblyIndex.load(
                    std::memory_order_relaxed
                );
                while (
                    assemblyIndex < previous &&
                    !parallelControl_->bestAssemblyIndex.compare_exchange_weak(
                        previous,
                        assemblyIndex,
                        std::memory_order_relaxed
                    )
                ) {}
            }
        }
        if (
            options_.targetAssemblyIndex >= 0 &&
            assemblyIndex <= options_.targetAssemblyIndex
        ) targetReached_ = true;
    }

    [[nodiscard]] std::u32string canonicalText(const Interval &interval) const
    {
        std::u32string text = original_.substr(
            static_cast<size_t>(interval.offset),
            static_cast<size_t>(interval.length)
        );
        if (!options_.acceptReversed) return text;

        std::u32string reversed(text.rbegin(), text.rend());
        return std::min(text, reversed);
    }

    int canonise(const Interval &interval, bool addToRollingHash = false)
    {
        auto intervalEntry = intervalIds_.find(interval);
        int id = -1;
        if (intervalEntry != intervalIds_.end()) id = intervalEntry->second;
        else
        {
            const std::u32string text = canonicalText(interval);
            auto textEntry = stringIds_.find(text);
            if (textEntry == stringIds_.end())
            {
                id = static_cast<int>(stringIds_.size());
                stringIds_.emplace(text, id);
            }
            else id = textEntry->second;
            intervalIds_.emplace(interval, id);
        }

        if (addToRollingHash)
            rollingHash_[id].push_back(interval.offset);
        return id;
    }

    [[nodiscard]] int knownCanonicalId(const Interval &interval) const
    {
        const auto entry = intervalIds_.find(interval);
        return entry == intervalIds_.end() ? -1 : entry->second;
    }

    [[nodiscard]] std::vector<Interval> preprocess() const
    {
        std::unordered_map<char32_t, int> firstOffsets;
        std::vector<bool> duplicated(original_.size(), false);
        for (size_t index = 0; index < original_.size(); index++)
        {
            const char32_t symbol = original_[index];
            const auto [entry, inserted] = firstOffsets.try_emplace(
                symbol,
                static_cast<int>(index)
            );
            if (!inserted)
            {
                duplicated[index] = true;
                duplicated[static_cast<size_t>(entry->second)] = true;
            }
        }

        FixedIntervalMap duplicateRuns(1);
        for (size_t index = 0; index < duplicated.size(); index++)
        {
            if (duplicated[index])
                duplicateRuns.insert(static_cast<int>(index));
        }
        return duplicateRuns.intervals(true);
    }

    [[nodiscard]] std::vector<int> stateKey(const AssemblyState &state)
    {
        std::vector<int> result;
        result.reserve(state.intervals.size());
        for (const Interval &interval : state.intervals)
            result.push_back(canonise(interval));
        if (result.size() > 1)
            std::sort(result.begin() + 1, result.end());
        return result;
    }

    Enumeration enumerate(AssemblyState &state, bool initial)
    {
        Enumeration result;
        const int ordinal = knownCanonicalId(state.intervals.front());
        const int maximumOrdinal = ordinal < 0
            ? std::numeric_limits<int>::max() : ordinal;
        std::vector<FixedIntervalMap> survivingIntervals(
            state.intervals.size()
        );

        std::vector<PotentialDuplicate> previous;
        for (size_t fragment = 0; fragment < state.intervals.size(); fragment++)
        {
            const Interval &interval = state.intervals[fragment];
            for (int offset = 0; offset < interval.length; offset++)
            {
                const Interval singleton{interval.offset + offset, 1};
                if (initial) canonise(singleton, true);
                PotentialDuplicate candidate{singleton, fragment};
                candidate.extend(previous, interval);
            }
        }

        bool active = true;
        bool overweight = false;
        size_t previousLength = 1;
        while (active)
        {
            if (shouldStop()) return {};
            result.duplicateSetsByLength.emplace_back();
            std::map<int, DuplicateSet> &sets =
                result.duplicateSetsByLength.back();
            active = false;
            std::vector<PotentialDuplicate> current;

            for (const PotentialDuplicate &candidate : previous)
            {
                if (shouldStop()) return {};
                const int id = canonise(candidate.interval, initial);
                if (id <= maximumOrdinal)
                {
                    auto [entry, inserted] = sets.try_emplace(
                        id,
                        previousLength + 1
                    );
                    static_cast<void>(inserted);
                    entry->second.insert(candidate);
                }
                else overweight = true;
            }

            for (const auto &[id, duplicates] : sets)
            {
                if (shouldStop()) return {};
                static_cast<void>(id);
                if (!duplicates.isValid()) continue;
                active = true;
                if (overweight || previousLength == 1)
                {
                    duplicates.extendValidOccurrences(
                        current,
                        state.intervals,
                        &survivingIntervals
                    );
                }
                else
                {
                    duplicates.extendValidOccurrences(
                        current,
                        state.intervals
                    );
                }
            }
            if (overweight) active = false;
            previousLength++;
            previous = std::move(current);
        }

        result.remnantIntervals.resize(survivingIntervals.size());
        for (size_t index = 0; index < survivingIntervals.size(); index++)
        {
            result.remnantIntervals[index] =
                survivingIntervals[index].intervals();
            for (const Interval &interval : result.remnantIntervals[index])
                canonise(interval);
        }
        if (initial)
        {
            for (auto &[id, offsets] : rollingHash_)
            {
                static_cast<void>(id);
                std::sort(offsets.begin(), offsets.end());
            }
        }
        return result;
    }

    [[nodiscard]] int lempelZivDuplicateBound(
        const AssemblyState &state
    ) const
    {
        std::vector<Interval> processedIntervals;
        int matches = 0;
        const int maximumFragmentLength = state.intervals.front().length;

        for (const Interval &currentInterval : state.intervals)
        {
            Interval window{currentInterval.offset, 1};
            int processedLimit = 0;
            for (int position = 1; position <= currentInterval.length; position++)
            {
                bool match = false;
                const int id = knownCanonicalId(window);
                const auto rollingEntry = rollingHash_.find(id);
                if (
                    id != -1 && rollingEntry != rollingHash_.end() &&
                    rollingEntry->second.size() > 1
                )
                {
                    const std::vector<int> &offsets = rollingEntry->second;
                    if (processedLimit > 0)
                    {
                        const auto occurrence = std::lower_bound(
                            offsets.begin(),
                            offsets.end(),
                            currentInterval.offset
                        );
                        if (
                            occurrence != offsets.end() &&
                            *occurrence + window.length <=
                                currentInterval.offset + processedLimit
                        ) match = true;
                    }
                    if (!match)
                    {
                        for (const Interval &processed : processedIntervals)
                        {
                            const auto occurrence = std::lower_bound(
                                offsets.begin(),
                                offsets.end(),
                                processed.offset
                            );
                            if (
                                occurrence != offsets.end() &&
                                *occurrence + window.length <=
                                    processed.offset + processed.length
                            )
                            {
                                match = true;
                                break;
                            }
                        }
                    }
                }

                if (match && window.length <= maximumFragmentLength)
                {
                    if (
                        window.offset + window.length ==
                        currentInterval.offset + currentInterval.length
                    )
                    {
                        matches += window.length - 1;
                        processedLimit = position;
                        window.offset = currentInterval.offset + position;
                        window.length = 1;
                    }
                    else window.length++;
                }
                else
                {
                    matches += window.length - 1;
                    if (window.length > 1)
                    {
                        matches--;
                        position--;
                    }
                    processedLimit = position;
                    window.offset = currentInterval.offset + position;
                    window.length = 1;
                }
            }
            processedIntervals.push_back(currentInterval);
        }
        return matches;
    }

    AssemblyState fragment(
        const ValidMatching &matching,
        const std::vector<std::vector<Interval>> &remnantIntervals
    )
    {
        Interval firstContainer;
        Interval secondContainer;
        bool foundFirst = false;
        bool foundSecond = false;
        for (const Interval &candidate : remnantIntervals[matching.firstFragment])
        {
            if (
                candidate.offset <= matching.first.offset &&
                candidate.offset + candidate.length >=
                    matching.first.offset + matching.first.length
            )
            {
                firstContainer = candidate;
                foundFirst = true;
                break;
            }
        }
        for (const Interval &candidate : remnantIntervals[matching.secondFragment])
        {
            if (
                candidate.offset <= matching.second.offset &&
                candidate.offset + candidate.length >=
                    matching.second.offset + matching.second.length
            )
            {
                secondContainer = candidate;
                foundSecond = true;
                break;
            }
        }
        if (!foundFirst || !foundSecond)
            throw std::logic_error("string matching escaped its remnant interval");

        AssemblyState result;
        result.intervals.push_back(matching.first);
        if (firstContainer == secondContainer)
        {
            Interval left = matching.first;
            Interval right = matching.second;
            if (right.offset < left.offset) std::swap(left, right);

            const int between = right.offset - (left.offset + left.length);
            if (between > 1)
                result.intervals.push_back(
                    {left.offset + left.length, between}
                );
            const int before = left.offset - firstContainer.offset;
            if (before > 1)
                result.intervals.push_back(
                    {firstContainer.offset, before}
                );
            const int after =
                firstContainer.offset + firstContainer.length -
                (right.offset + right.length);
            if (after > 1)
                result.intervals.push_back(
                    {right.offset + right.length, after}
                );
        }
        else
        {
            const int beforeFirst =
                matching.first.offset - firstContainer.offset;
            if (beforeFirst > 1)
                result.intervals.push_back(
                    {firstContainer.offset, beforeFirst}
                );
            const int afterFirst =
                firstContainer.offset + firstContainer.length -
                (matching.first.offset + matching.first.length);
            if (afterFirst > 1)
                result.intervals.push_back(
                    {
                        matching.first.offset + matching.first.length,
                        afterFirst
                    }
                );
            const int beforeSecond =
                matching.second.offset - secondContainer.offset;
            if (beforeSecond > 1)
                result.intervals.push_back(
                    {secondContainer.offset, beforeSecond}
                );
            const int afterSecond =
                secondContainer.offset + secondContainer.length -
                (matching.second.offset + matching.second.length);
            if (afterSecond > 1)
                result.intervals.push_back(
                    {
                        matching.second.offset + matching.second.length,
                        afterSecond
                    }
                );
        }

        for (const Interval &interval : result.intervals) canonise(interval);
        for (const std::vector<Interval> &fragmentRemnants : remnantIntervals)
        {
            for (const Interval &interval : fragmentRemnants)
            {
                if (
                    interval != firstContainer &&
                    interval != secondContainer
                ) result.intervals.push_back(interval);
            }
        }
        return result;
    }

    /** Explore one match using this worker's canonical IDs, cache and pathway. */
    void exploreMatching(
        const ValidMatching &matching,
        const std::vector<std::vector<Interval>> &remnantIntervals,
        int parentDuplicatedSymbols
    )
    {
        AssemblyState next = fragment(matching, remnantIntervals);
        next.duplicatedSymbols =
            parentDuplicatedSymbols + matching.fragmentLength - 1;

        const int lowerBound =
            static_cast<int>(original_.size()) -
            next.duplicatedSymbols - 1 - lempelZivDuplicateBound(next);
        if (lowerBound >= pruningIndex()) return;

        std::vector<int> key = stateKey(next);
        const auto existing = states_.find(key);
        if (
            existing != states_.end() &&
            existing->second >= next.duplicatedSymbols
        ) return;

        states_.insert_or_assign(std::move(key), next.duplicatedSymbols);
        currentPath_.push_back({matching.first, matching.second});
        recurse(next, false);
        currentPath_.pop_back();
    }

    void recurse(AssemblyState &state, bool initial)
    {
        if (shouldStop()) return;

        const int assemblyIndex =
            static_cast<int>(original_.size()) - state.duplicatedSymbols - 1;
        recordBest(assemblyIndex, currentPath_);
        if (targetReached_ || optimumReached()) return;

        Enumeration enumeration = enumerate(state, initial);
        if (shouldStop()) return;
        for (
            size_t levelIndex = enumeration.duplicateSetsByLength.size();
            levelIndex-- > 0;
        )
        {
            std::map<int, DuplicateSet> &sets =
                enumeration.duplicateSetsByLength[levelIndex];
            for (const auto &[id, duplicateSet] : sets)
            {
                static_cast<void>(id);
                auto matchings = duplicateSet.matchingCursor();
                ValidMatching matching;
                while (matchings.next(matching))
                {
                    if (shouldStop()) return;
                    exploreMatching(
                        matching,
                        enumeration.remnantIntervals,
                        state.duplicatedSymbols
                    );
                }
            }
        }
    }

    void runParallel(AssemblyState &root)
    {
        // Canonical ordinals constrain descendant enumeration. Establish them
        // once in serial order, then give every worker the same seed. Sharing
        // a mutable ID allocator would make pruning depend on scheduling.
        const Enumeration enumeration = enumerate(root, true);
        if (shouldStop()) return;

        // Keep the serial root order and shard ordinals without storing every
        // pair. Only cursor advancement is shared; descendants remain local.
        size_t levelIndex = enumeration.duplicateSetsByLength.size();
        const std::map<int, DuplicateSet> *sets = nullptr;
        std::map<int, DuplicateSet>::const_iterator nextSet;
        DuplicateSet::MatchingCursor matchingCursor;
        size_t ordinal = 0;
        const auto nextMatching = [&](Search &search, ValidMatching &matching)
        {
            for (;;)
            {
                if (search.shouldStop()) return false;
                if (matchingCursor.next(matching))
                {
                    const bool selected =
                        ordinal % options_.shardCount == options_.shardIndex;
                    ordinal++;
                    if (selected) return true;
                    continue;
                }

                while (sets == nullptr || nextSet == sets->end())
                {
                    if (levelIndex == 0) return false;
                    sets = &enumeration.duplicateSetsByLength[--levelIndex];
                    nextSet = sets->begin();
                }
                matchingCursor = (nextSet++)->second.matchingCursor();
            }
        };

        // Prefetch at most one job per requested worker so a small search does
        // not start unnecessary threads. All workers share these jobs because
        // the OpenMP runtime may supply fewer threads than requested.
        std::vector<ValidMatching> initialJobs;
        for (int workerIndex = 0; workerIndex < options_.threadCount; workerIndex++)
        {
            ValidMatching matching;
            if (!nextMatching(*this, matching)) break;
            initialJobs.push_back(matching);
        }
        if (initialJobs.empty() || shouldStop()) return;

        ParallelControl control(bestAssemblyIndex_);
        std::mutex jobsMutex;
        size_t nextInitialJob = 0;
        const int threadCount = static_cast<int>(initialJobs.size());
        const auto acquireJob = [&](Search &worker, ValidMatching &matching)
        {
            std::lock_guard<std::mutex> lock(jobsMutex);
            if (nextInitialJob < initialJobs.size())
            {
                matching = initialJobs[nextInitialJob++];
                return true;
            }
            return nextMatching(worker, matching);
        };
        std::vector<Result> results(static_cast<size_t>(threadCount));
        const auto work = [&](int workerIndex)
        {
            try
            {
                Search worker = *this;
                worker.parallelControl_ = &control;
                for (;;)
                {
                    if (worker.shouldStop()) break;
                    ValidMatching matching;
                    if (!acquireJob(worker, matching)) break;
                    worker.exploreMatching(matching, enumeration.remnantIntervals, 0);
                }
                Result &result = results[static_cast<size_t>(workerIndex)];
                result.assemblyIndex = worker.bestAssemblyIndex_;
                result.pathway = std::move(worker.bestPath_);
            }
            catch (...)
            {
                std::lock_guard<std::mutex> lock(control.failureMutex);
                if (control.failure == nullptr) control.failure = std::current_exception();
                control.failed.store(true, std::memory_order_relaxed);
            }
        };

#if defined(PARALLELASSEMBLYCPP_USE_OPENMP)
        if (threadCount > 1)
        {
            #pragma omp parallel num_threads(threadCount)
            {
                work(omp_get_thread_num());
            }
        }
        else
#endif
            work(0);

        if (control.failure != nullptr) std::rethrow_exception(control.failure);
        runtimeLimitReached_ = control.runtimeLimitReached.load(std::memory_order_relaxed);
        interrupted_ = control.interrupted.load(std::memory_order_relaxed);
        for (Result &result : results)
        {
            // A runtime may provide fewer threads than requested. Unused
            // result slots have the default sentinel and no witness.
            if (result.assemblyIndex < 0) continue;
            if (result.assemblyIndex < bestAssemblyIndex_)
            {
                bestAssemblyIndex_ = result.assemblyIndex;
                bestPath_ = std::move(result.pathway);
            }
        }
    }

    void runFromRoot()
    {
        AssemblyState root;
        if (shouldStop()) return;
        recordBest(bestAssemblyIndex_, bestPath_);
        if (targetReached_) return;
        root.intervals = preprocess();
        if (shouldStop() || root.intervals.empty()) return;
        try
        {
            const auto stop = [this] { return shouldStop(); };
            const auto grammar = stringRepair::implementation::Compressor(
                original_, options_.acceptReversed, stop
            ).run();
            const auto pathway = repairPathway(grammar, stop);
            int index = static_cast<int>(original_.size()) - 1;
            for (const auto &step : pathway) index -= step.match.length - 1;
            // The path must accompany the bound: >= pruning can discard every
            // exact-search branch tying an already optimal Re-Pair incumbent.
            recordBest(index, pathway);
        }
        catch (const stringRepair::implementation::Cancelled &)
        {
            return;
        }
        if (shouldStop()) return;
        if (assemblyLowerBound_ < 0)
        {
            assemblyLowerBound_ = assembly_bounds::scalarLowerBound(
                static_cast<int>(original_.size())
            );
            if (assembly_bounds::vectorBoundsEnabled && !optimumReached())
            {
                // Concatenation adds symbol counts, including when either
                // operand is reversed. Every string construction therefore
                // maps to a vector addition chain. Its certified lower bound
                // can prove a complete witness already held here optimal.
                std::unordered_map<char32_t, int> frequencies;
                for (size_t index = 0; index < original_.size(); ++index)
                {
                    if ((index & 1023U) == 0 && shouldStop()) return;
                    ++frequencies[original_[index]];
                }
                std::vector<int> counts;
                counts.reserve(frequencies.size());
                for (const auto &[symbol, count] : frequencies)
                {
                    static_cast<void>(symbol);
                    counts.push_back(count);
                }
                assembly_bounds::VectorBoundCache bounds;
                assemblyLowerBound_ = bounds.lowerBound(
                    std::move(counts), bestAssemblyIndex_
                );
            }
        }
        if (shouldStop()) return;
        std::vector<int> rootKey(root.intervals.size(), -1);
        states_.emplace(std::move(rootKey), 0);
        if (options_.threadCount > 1 || options_.shardCount > 1)
            runParallel(root);
        else
            recurse(root, true);
    }

public:
    Search(std::u32string original, const Options &options):
        original_(std::move(original)), options_(options) {}

    Result run()
    {
        if (
            original_.size() >=
            static_cast<size_t>(std::numeric_limits<int>::max())
        ) throw std::invalid_argument("string is too long to index");
        if (options_.threadCount < 1)
            throw std::invalid_argument("string thread count must be positive");
        if (options_.shardCount == 0 || options_.shardIndex >= options_.shardCount)
            throw std::invalid_argument("invalid string search shard");
#if !defined(PARALLELASSEMBLYCPP_USE_OPENMP)
        if (options_.threadCount > 1)
            throw std::invalid_argument("parallel string threads require an OpenMP-enabled executable");
#endif

        started_ = std::clock();
        bestAssemblyIndex_ = static_cast<int>(original_.size()) - 1;
        runFromRoot();

        if (
            (options_.threadCount > 1 || options_.shardCount > 1) &&
            options_.reconstructPathway && !runtimeLimitReached_ && !interrupted_ &&
            bestAssemblyIndex_ >= 0 && !bestPath_.empty()
        )
        {
            Options replayOptions = options_;
            replayOptions.threadCount = 1;
            replayOptions.shardIndex = 0;
            replayOptions.shardCount = 1;
            replayOptions.targetAssemblyIndex = bestAssemblyIndex_;
            replayOptions.reconstructPathway = false;
            Search replay(original_, replayOptions);
            replay.started_ = started_;
            replay.bestAssemblyIndex_ = static_cast<int>(original_.size()) - 1;
            replay.assemblyLowerBound_ = assemblyLowerBound_;
            replay.runFromRoot();
            if (replay.bestAssemblyIndex_ <= bestAssemblyIndex_)
            {
                bestAssemblyIndex_ = replay.bestAssemblyIndex_;
                bestPath_ = std::move(replay.bestPath_);
            }
            runtimeLimitReached_ = replay.runtimeLimitReached_;
            interrupted_ = replay.interrupted_;
        }

        Result result;
        result.assemblyIndex = bestAssemblyIndex_;
        result.clockTicks = elapsedTicks();
        result.runtimeLimitReached = runtimeLimitReached_;
        result.interrupted = interrupted_;
        result.pathway = bestPath_;
        return result;
    }
};

inline std::vector<Interval> remnantIntervals(
    size_t stringLength,
    const std::vector<PathwayStep> &pathway
)
{
    std::vector<Interval> removed;
    removed.reserve(pathway.size());
    for (const PathwayStep &step : pathway)
        removed.push_back(step.duplicate);
    std::sort(
        removed.begin(),
        removed.end(),
        [](const Interval &left, const Interval &right)
        {
            return left.offset < right.offset;
        }
    );

    std::vector<Interval> merged;
    for (const Interval &interval : removed)
    {
        if (
            merged.empty() ||
            merged.back().offset + merged.back().length < interval.offset
        ) merged.push_back(interval);
        else
        {
            const int end = std::max(
                merged.back().offset + merged.back().length,
                interval.offset + interval.length
            );
            merged.back().length = end - merged.back().offset;
        }
    }

    std::vector<Interval> result;
    int cursor = 0;
    for (const Interval &interval : merged)
    {
        if (cursor < interval.offset)
            result.push_back({cursor, interval.offset - cursor});
        cursor = std::max(cursor, interval.offset + interval.length);
    }
    const int length = static_cast<int>(stringLength);
    if (cursor < length) result.push_back({cursor, length - cursor});
    return result;
}

} // namespace implementation

/**
 * Search one UTF-8 string, counting Unicode scalar values without normalization.
 * A stopped or sharded search returns its best witnessed index, which need not
 * be the global minimum. clockTicks excludes UTF-8 decoding.
 *
 * @throws std::invalid_argument for invalid UTF-8 or unsupported search options
 */
inline Result calculate(const std::string &input, const Options &options = {})
{
    return implementation::Search(implementation::decodeInput(input), options).run();
}

/** Write ASCII JSON for a witness belonging to input; report failures in error. */
inline bool writePathway(
    const std::string &filename,
    const std::string &input,
    const Result &result,
    std::string &error
)
{
    std::ofstream output;
    try
    {
        std::vector<size_t> byteOffsets;
        implementation::decodeInput(input, &byteOffsets);
        output.open(filename);
        if (!output.is_open())
        {
            error = "could not open output file '" + filename + "'";
            return false;
        }
        output << "{\n  \"file_graph\": [\n    {\n      \"Fragments\": [";
        implementation::writeJsonString(input, output);
        output << "],\n      \"Positions\": [0]\n    }\n  ],\n";

        const std::vector<Interval> remnants =
            implementation::remnantIntervals(byteOffsets.size() - 1, result.pathway);
        output << "  \"remnant\": [\n    {\n      \"Fragments\": [";
        for (size_t index = 0; index < remnants.size(); index++)
        {
            if (index > 0) output << ',';
            const size_t first = byteOffsets.at(static_cast<size_t>(remnants[index].offset));
            const size_t last = byteOffsets.at(static_cast<size_t>(
                remnants[index].offset + remnants[index].length
            ));
            implementation::writeJsonString(
                std::string_view(input).substr(first, last - first), output
            );
        }
        output << "],\n      \"Positions\": [";
        for (size_t index = 0; index < remnants.size(); index++)
        {
            if (index > 0) output << ',';
            output << remnants[index].offset;
        }
        output << "]\n    }\n  ],\n  \"duplicates\": [";
        for (size_t index = 0; index < result.pathway.size(); index++)
        {
            if (index > 0) output << ',';
            const PathwayStep &step = result.pathway[index];
            output << "\n    {\"Left\":[" << step.match.offset << ','
                   << step.match.length << "],\"Right\":["
                   << step.duplicate.offset << ',' << step.duplicate.length
                   << "]}";
        }
        if (!result.pathway.empty()) output << '\n';
        output << "  ]\n}\n";
    }
    catch (const std::exception &failure)
    {
        error = "could not write output file '" + filename + "': " +
            failure.what();
        return false;
    }
    output.close();
    if (!output)
    {
        error = "could not write output file '" + filename + "'";
        return false;
    }
    return true;
}

} // namespace parallelassemblycpp::detail::stringAssembly
