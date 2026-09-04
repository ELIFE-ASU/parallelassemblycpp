#pragma once

/*
 * Adapted for AssemblyCpp v5 from the string assembly implementation in
 * croningroup/public/assemblycpp-public commit 2a87948, authored by Stuart
 * Marshall from work by Ian Seet and Leroy Cronin. See README.md and
 * License.md for source details and licensing.
 */

#include <algorithm>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace assemblycpp::detail::stringAssembly
{

/** A half-open interval in the original string. */
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
};

struct Result
{
    int assemblyIndex = -1;
    unsigned long long clockTicks = 0;
    bool runtimeLimitReached = false;
    bool interrupted = false;
    std::vector<PathwayStep> pathway;
};

namespace implementation
{

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
    std::vector<FixedIntervalMap> intervalsByFragment_;
    std::vector<PotentialDuplicate> occurrences_;

public:
    DuplicateSet(size_t fragmentLength, size_t fragmentCount):
        fragmentLength_(fragmentLength),
        intervalsByFragment_(
            fragmentCount,
            FixedIntervalMap(static_cast<int>(fragmentLength))
        ) {}

    void insert(const PotentialDuplicate &occurrence)
    {
        occurrences_.push_back(occurrence);
        intervalsByFragment_[occurrence.fragment].insert(
            occurrence.interval.offset
        );
    }

    [[nodiscard]] bool isValid() const
    {
        size_t populatedFragments = 0;
        const FixedIntervalMap *onlyFragment = nullptr;
        for (const FixedIntervalMap &intervals : intervalsByFragment_)
        {
            if (intervals.empty()) continue;
            populatedFragments++;
            onlyFragment = &intervals;
            if (populatedFragments > 1) return true;
        }
        return onlyFragment != nullptr &&
            onlyFragment->containsTwoDisjointInsertions();
    }

    [[nodiscard]] std::vector<ValidMatching> matchings() const
    {
        std::vector<ValidMatching> result;
        for (size_t first = 0; first < occurrences_.size(); first++)
        {
            for (
                size_t second = first + 1;
                second < occurrences_.size();
                second++
            )
            {
                const PotentialDuplicate &left = occurrences_[first];
                const PotentialDuplicate &right = occurrences_[second];
                if (
                    left.fragment == right.fragment &&
                    left.overlaps(right)
                ) continue;
                result.push_back(
                    {
                        left.interval,
                        right.interval,
                        left.fragment,
                        right.fragment,
                        static_cast<int>(fragmentLength_)
                    }
                );
            }
        }
        return result;
    }

    bool extendValidOccurrences(
        std::vector<PotentialDuplicate> &output,
        const std::vector<Interval> &fragments,
        std::vector<FixedIntervalMap> *survivingIntervals = nullptr
    ) const
    {
        std::vector<bool> valid(occurrences_.size(), false);
        for (size_t first = 0; first < occurrences_.size(); first++)
        {
            for (
                size_t second = first + 1;
                second < occurrences_.size();
                second++
            )
            {
                const PotentialDuplicate &left = occurrences_[first];
                const PotentialDuplicate &right = occurrences_[second];
                if (
                    left.fragment == right.fragment &&
                    left.overlaps(right)
                ) continue;

                valid[first] = true;
                valid[second] = true;
                if (survivingIntervals != nullptr)
                {
                    (*survivingIntervals)[left.fragment].insert(
                        left.interval.offset
                    );
                    (*survivingIntervals)[right.fragment].insert(
                        right.interval.offset
                    );
                }
            }
        }

        bool extendedAny = false;
        for (size_t index = 0; index < occurrences_.size(); index++)
        {
            if (!valid[index]) continue;
            occurrences_[index].extend(
                output,
                fragments[occurrences_[index].fragment]
            );
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

class Search
{
    std::string original_;
    Options options_;
    std::clock_t started_ = 0;
    std::unordered_map<std::string, int> stringIds_;
    std::unordered_map<Interval, int, IntervalHash> intervalIds_;
    std::unordered_map<int, std::vector<int>> rollingHash_;
    std::unordered_map<std::vector<int>, int, IntegerVectorHash> states_;
    std::vector<PathwayStep> currentPath_;
    std::vector<PathwayStep> bestPath_;
    int bestAssemblyIndex_ = -1;
    bool runtimeLimitReached_ = false;
    bool interrupted_ = false;

    [[nodiscard]] unsigned long long elapsedTicks() const noexcept
    {
        const std::clock_t now = std::clock();
        const std::clock_t error = static_cast<std::clock_t>(-1);
        if (now == error || started_ == error) return 0;

        using UnsignedClock = std::make_unsigned_t<std::clock_t>;
        const UnsignedClock elapsed =
            static_cast<UnsignedClock>(now) -
            static_cast<UnsignedClock>(started_);
        if constexpr (sizeof(UnsignedClock) > sizeof(unsigned long long))
        {
            if (
                elapsed >
                static_cast<UnsignedClock>(
                    std::numeric_limits<unsigned long long>::max()
                )
            ) return std::numeric_limits<unsigned long long>::max();
        }
        return static_cast<unsigned long long>(elapsed);
    }

    bool shouldStop()
    {
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

    [[nodiscard]] std::string canonicalText(const Interval &interval) const
    {
        std::string text = original_.substr(
            static_cast<size_t>(interval.offset),
            static_cast<size_t>(interval.length)
        );
        if (!options_.acceptReversed) return text;

        std::string reversed(text.rbegin(), text.rend());
        return std::min(text, reversed);
    }

    int canonise(const Interval &interval, bool addToRollingHash = false)
    {
        auto intervalEntry = intervalIds_.find(interval);
        int id = -1;
        if (intervalEntry != intervalIds_.end()) id = intervalEntry->second;
        else
        {
            const std::string text = canonicalText(interval);
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
        std::unordered_map<unsigned char, int> firstOffsets;
        std::vector<bool> duplicated(original_.size(), false);
        for (size_t index = 0; index < original_.size(); index++)
        {
            const unsigned char symbol =
                static_cast<unsigned char>(original_[index]);
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
            result.duplicateSetsByLength.emplace_back();
            std::map<int, DuplicateSet> &sets =
                result.duplicateSetsByLength.back();
            active = false;
            std::vector<PotentialDuplicate> current;

            for (const PotentialDuplicate &candidate : previous)
            {
                const int id = canonise(candidate.interval, initial);
                if (id <= maximumOrdinal)
                {
                    auto [entry, inserted] = sets.try_emplace(
                        id,
                        previousLength + 1,
                        state.intervals.size()
                    );
                    static_cast<void>(inserted);
                    entry->second.insert(candidate);
                }
                else overweight = true;
            }

            for (const auto &[id, duplicates] : sets)
            {
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

    void recurse(AssemblyState &state, bool initial)
    {
        if (shouldStop()) return;

        const int assemblyIndex =
            static_cast<int>(original_.size()) - state.duplicatedSymbols - 1;
        if (assemblyIndex < bestAssemblyIndex_)
        {
            bestAssemblyIndex_ = assemblyIndex;
            bestPath_ = currentPath_;
        }

        Enumeration enumeration = enumerate(state, initial);
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
                std::vector<ValidMatching> matchings =
                    duplicateSet.matchings();
                for (size_t matchingIndex = matchings.size(); matchingIndex-- > 0;)
                {
                    if (shouldStop()) return;
                    const ValidMatching &matching = matchings[matchingIndex];
                    AssemblyState next = fragment(
                        matching,
                        enumeration.remnantIntervals
                    );
                    next.duplicatedSymbols =
                        state.duplicatedSymbols + matching.fragmentLength - 1;

                    const int lowerBound =
                        static_cast<int>(original_.size()) -
                        next.duplicatedSymbols - 1 -
                        lempelZivDuplicateBound(next);
                    if (lowerBound >= bestAssemblyIndex_) continue;

                    std::vector<int> key = stateKey(next);
                    const PathwayStep step{matching.first, matching.second};
                    const auto existing = states_.find(key);
                    if (
                        existing == states_.end() ||
                        next.duplicatedSymbols > existing->second
                    )
                    {
                        states_.insert_or_assign(
                            std::move(key),
                            next.duplicatedSymbols
                        );
                        currentPath_.push_back(step);
                        recurse(next, false);
                        currentPath_.pop_back();
                    }
                }
            }
        }
    }

public:
    Search(std::string original, const Options &options):
        original_(std::move(original)), options_(options) {}

    Result run()
    {
        if (
            original_.size() >=
            static_cast<size_t>(std::numeric_limits<int>::max())
        ) throw std::invalid_argument("string is too long to index");

        started_ = std::clock();
        bestAssemblyIndex_ = static_cast<int>(original_.size()) - 1;

        AssemblyState root;
        root.intervals = preprocess();
        std::vector<int> rootKey(root.intervals.size(), -1);
        states_.emplace(std::move(rootKey), 0);

        if (!root.intervals.empty()) recurse(root, true);

        Result result;
        result.assemblyIndex = bestAssemblyIndex_;
        result.clockTicks = elapsedTicks();
        result.runtimeLimitReached = runtimeLimitReached_;
        result.interrupted = interrupted_;
        result.pathway = bestPath_;
        return result;
    }
};

inline void writeJsonString(std::string_view value, std::ostream &output)
{
    static constexpr char hexDigits[] = "0123456789ABCDEF";
    output.put('"');
    for (const unsigned char character : value)
    {
        switch (character)
        {
            case '"': output << "\\\""; break;
            case '\\': output << "\\\\"; break;
            case '\b': output << "\\b"; break;
            case '\f': output << "\\f"; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default:
                if (character < 0x20)
                {
                    output << "\\u00"
                           << hexDigits[character >> 4]
                           << hexDigits[character & 0x0f];
                }
                else output.put(static_cast<char>(character));
                break;
        }
    }
    output.put('"');
}

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

inline Result calculate(const std::string &input, const Options &options = {})
{
    return implementation::Search(input, options).run();
}

inline bool writePathway(
    const std::string &filename,
    const std::string &input,
    const Result &result,
    std::string &error
)
{
    std::ofstream output(filename);
    if (!output.is_open())
    {
        error = "could not open output file '" + filename + "'";
        return false;
    }

    output << "{\n  \"file_graph\": [\n    {\n      \"Fragments\": [";
    implementation::writeJsonString(input, output);
    output << "],\n      \"Positions\": [0]\n    }\n  ],\n";

    const std::vector<Interval> remnants =
        implementation::remnantIntervals(input.size(), result.pathway);
    output << "  \"remnant\": [\n    {\n      \"Fragments\": [";
    for (size_t index = 0; index < remnants.size(); index++)
    {
        if (index > 0) output << ',';
        implementation::writeJsonString(
            std::string_view(input).substr(
                static_cast<size_t>(remnants[index].offset),
                static_cast<size_t>(remnants[index].length)
            ),
            output
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
    output.close();
    if (!output)
    {
        error = "could not write output file '" + filename + "'";
        return false;
    }
    return true;
}

} // namespace assemblycpp::detail::stringAssembly
