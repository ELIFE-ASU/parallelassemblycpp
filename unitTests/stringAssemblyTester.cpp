#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../v5/stringAssembly.h"

using assemblycpp::detail::stringAssembly::Interval;
using assemblycpp::detail::stringAssembly::Options;
using assemblycpp::detail::stringAssembly::PathwayStep;
using assemblycpp::detail::stringAssembly::Result;
using assemblycpp::detail::stringAssembly::calculate;
using assemblycpp::detail::stringAssembly::writePathway;

namespace implementation =
    assemblycpp::detail::stringAssembly::implementation;

namespace
{

class TestFailure: public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

void require(bool condition, const std::string &message)
{
    if (!condition) throw TestFailure(message);
}

std::string repeated(char symbol, size_t count)
{
    return std::string(count, symbol);
}

std::string numberedBlocks(size_t blockLength)
{
    return repeated('0', blockLength) +
        repeated('1', blockLength) +
        repeated('2', blockLength);
}

std::string intervalListText(const std::vector<Interval> &intervals)
{
    std::ostringstream output;
    output << '[';
    for (size_t index = 0; index < intervals.size(); index++)
    {
        if (index > 0) output << ", ";
        output << '{' << intervals[index].offset << ','
               << intervals[index].length << '}';
    }
    output << ']';
    return output.str();
}

std::string withoutJsonFormatting(const std::string &json)
{
    std::string result;
    result.reserve(json.size());
    bool insideString = false;
    bool escaped = false;
    for (const char character : json)
    {
        if (
            !insideString &&
            (character == ' ' || character == '\t' || character == '\n' ||
                character == '\r')
        ) continue;

        result.push_back(character);
        if (!insideString)
        {
            if (character == '"') insideString = true;
            continue;
        }
        if (escaped) escaped = false;
        else if (character == '\\') escaped = true;
        else if (character == '"') insideString = false;
    }
    return result;
}

void requireIntervals(
    const std::vector<Interval> &actual,
    const std::vector<Interval> &expected,
    const std::string &description
)
{
    require(
        actual == expected,
        description + ": expected " + intervalListText(expected) +
            ", got " + intervalListText(actual)
    );
}

/**
 * Small, deliberately straightforward reference search. It models a state as
 * a multiset of fragments and tries every pair of equal, non-overlapping
 * substrings. This is independent of the production search's enumeration,
 * hashing, bounds, and pruning, making it useful as an exhaustive oracle for
 * short inputs.
 */
class ReferenceSearch
{
    struct Occurrence
    {
        size_t fragment = 0;
        size_t offset = 0;
        size_t length = 0;
    };

    bool acceptReversed_ = false;
    std::unordered_map<std::string, int> memo_;

    [[nodiscard]] std::string canonicalFragment(std::string fragment) const
    {
        if (!acceptReversed_) return fragment;
        std::string reversed(fragment.rbegin(), fragment.rend());
        return std::min(fragment, reversed);
    }

    [[nodiscard]] std::string stateKey(
        const std::vector<std::string> &fragments
    ) const
    {
        std::vector<std::string> canonical;
        canonical.reserve(fragments.size());
        for (const std::string &fragment : fragments)
            canonical.push_back(canonicalFragment(fragment));
        std::sort(canonical.begin(), canonical.end());

        std::string key;
        for (const std::string &fragment : canonical)
        {
            key += std::to_string(fragment.size());
            key.push_back(':');
            key += fragment;
            key.push_back(';');
        }
        return key;
    }

    [[nodiscard]] bool equivalent(
        const std::vector<std::string> &fragments,
        const Occurrence &left,
        const Occurrence &right
    ) const
    {
        const std::string leftText = fragments[left.fragment].substr(
            left.offset,
            left.length
        );
        const std::string rightText = fragments[right.fragment].substr(
            right.offset,
            right.length
        );
        if (leftText == rightText) return true;
        return acceptReversed_ && std::equal(
            leftText.begin(),
            leftText.end(),
            rightText.rbegin()
        );
    }

    static void appendPart(
        std::vector<std::string> &output,
        const std::string &fragment,
        size_t offset,
        size_t length
    )
    {
        if (length > 1) output.push_back(fragment.substr(offset, length));
    }

    [[nodiscard]] static std::vector<std::string> fragmentState(
        const std::vector<std::string> &fragments,
        const Occurrence &left,
        const Occurrence &right
    )
    {
        std::vector<std::string> result{
            fragments[left.fragment].substr(left.offset, left.length)
        };

        for (size_t index = 0; index < fragments.size(); index++)
        {
            const std::string &fragment = fragments[index];
            if (index != left.fragment && index != right.fragment)
            {
                if (fragment.size() > 1) result.push_back(fragment);
                continue;
            }

            if (left.fragment == right.fragment)
            {
                if (index != left.fragment) continue;
                const size_t firstOffset = std::min(left.offset, right.offset);
                const size_t secondOffset = std::max(left.offset, right.offset);
                appendPart(result, fragment, 0, firstOffset);
                const size_t betweenOffset = firstOffset + left.length;
                appendPart(
                    result,
                    fragment,
                    betweenOffset,
                    secondOffset - betweenOffset
                );
                const size_t afterOffset = secondOffset + left.length;
                appendPart(
                    result,
                    fragment,
                    afterOffset,
                    fragment.size() - afterOffset
                );
                continue;
            }

            const Occurrence &selected =
                index == left.fragment ? left : right;
            appendPart(result, fragment, 0, selected.offset);
            const size_t afterOffset = selected.offset + selected.length;
            appendPart(
                result,
                fragment,
                afterOffset,
                fragment.size() - afterOffset
            );
        }
        return result;
    }

    [[nodiscard]] int maximumAdditionalSavings(
        const std::vector<std::string> &fragments
    )
    {
        const std::string key = stateKey(fragments);
        const auto known = memo_.find(key);
        if (known != memo_.end()) return known->second;

        size_t maximumLength = 0;
        for (const std::string &fragment : fragments)
            maximumLength = std::max(maximumLength, fragment.size());

        int best = 0;
        std::unordered_set<std::string> triedTransitions;
        for (size_t length = 2; length <= maximumLength; length++)
        {
            std::vector<Occurrence> occurrences;
            for (size_t fragment = 0; fragment < fragments.size(); fragment++)
            {
                const size_t fragmentLength = fragments[fragment].size();
                if (fragmentLength < length) continue;
                for (size_t offset = 0; offset <= fragmentLength - length; offset++)
                    occurrences.push_back({fragment, offset, length});
            }

            for (size_t first = 0; first < occurrences.size(); first++)
            {
                for (size_t second = first + 1; second < occurrences.size(); second++)
                {
                    const Occurrence &left = occurrences[first];
                    const Occurrence &right = occurrences[second];
                    if (
                        left.fragment == right.fragment &&
                        left.offset < right.offset + length &&
                        right.offset < left.offset + length
                    ) continue;
                    if (!equivalent(fragments, left, right)) continue;

                    std::vector<std::string> next = fragmentState(
                        fragments,
                        left,
                        right
                    );
                    const std::string transitionKey =
                        std::to_string(length) + '#' + stateKey(next);
                    if (!triedTransitions.insert(transitionKey).second) continue;
                    best = std::max(
                        best,
                        static_cast<int>(length) - 1 +
                            maximumAdditionalSavings(next)
                    );
                }
            }
        }

        memo_.emplace(key, best);
        return best;
    }

public:
    explicit ReferenceSearch(bool acceptReversed):
        acceptReversed_(acceptReversed) {}

    [[nodiscard]] int assemblyIndex(const std::string &input)
    {
        return static_cast<int>(input.size()) - 1 -
            maximumAdditionalSavings({input});
    }
};

std::string intervalText(const std::string &input, const Interval &interval)
{
    return input.substr(
        static_cast<size_t>(interval.offset),
        static_cast<size_t>(interval.length)
    );
}

void requireValidInterval(
    const Interval &interval,
    size_t inputLength,
    const std::string &description
)
{
    require(interval.offset >= 0, description + " has a negative offset");
    require(interval.length > 0, description + " is empty");
    const size_t offset = static_cast<size_t>(interval.offset);
    const size_t length = static_cast<size_t>(interval.length);
    require(offset <= inputLength, description + " starts past the input");
    require(
        length <= inputLength - offset,
        description + " extends past the input"
    );
}

bool overlap(const Interval &left, const Interval &right)
{
    return left.offset < right.offset + right.length &&
        right.offset < left.offset + left.length;
}

void requireConsistentPathway(
    const std::string &input,
    const Result &result,
    bool acceptReversed = false
)
{
    int duplicatedSymbols = 0;
    std::vector<Interval> removed;
    for (size_t index = 0; index < result.pathway.size(); index++)
    {
        const PathwayStep &step = result.pathway[index];
        const std::string prefix =
            "pathway step " + std::to_string(index) + ' ';
        requireValidInterval(step.match, input.size(), prefix + "match");
        requireValidInterval(step.duplicate, input.size(), prefix + "duplicate");
        require(
            step.match.length == step.duplicate.length,
            prefix + "uses unequal fragment lengths"
        );
        require(step.match.length >= 2, prefix + "does not save an operation");
        require(!overlap(step.match, step.duplicate), prefix + "overlaps itself");
        for (const Interval &previous : removed)
        {
            require(
                !overlap(step.match, previous),
                prefix + "match reuses an already removed region"
            );
            require(
                !overlap(step.duplicate, previous),
                prefix + "duplicate reuses an already removed region"
            );
        }

        const std::string match = intervalText(input, step.match);
        const std::string duplicate = intervalText(input, step.duplicate);
        std::string reversedDuplicate(duplicate.rbegin(), duplicate.rend());
        require(
            match == duplicate || (acceptReversed && match == reversedDuplicate),
            prefix + "does not identify equivalent text"
        );
        duplicatedSymbols += step.match.length - 1;
        removed.push_back(step.duplicate);
    }

    require(
        result.assemblyIndex ==
            static_cast<int>(input.size()) - duplicatedSymbols - 1,
        "assembly index is inconsistent with the returned pathway"
    );
}

void testUpstreamRegressionCases()
{
    struct TestCase
    {
        size_t blockLength;
        int expectedIndex;
    };

    // These are the five strings and reference results from
    // assemblycpp-public/tests/strings/012test(_standardOutput).
    const std::vector<TestCase> cases{
        {5, 11},
        {10, 14},
        {15, 17},
        {20, 17},
        {25, 20}
    };

    for (const TestCase &testCase : cases)
    {
        const std::string input = numberedBlocks(testCase.blockLength);
        const Result result = calculate(input);
        require(
            result.assemblyIndex == testCase.expectedIndex,
            "upstream block-length " + std::to_string(testCase.blockLength) +
                " case: expected index " +
                std::to_string(testCase.expectedIndex) + ", got " +
                std::to_string(result.assemblyIndex)
        );
        require(!result.runtimeLimitReached, "unlimited search timed out");
        require(!result.interrupted, "unlimited search was interrupted");
        requireConsistentPathway(input, result);
    }
}

void testTrivialStrings()
{
    const Result empty = calculate("");
    require(empty.assemblyIndex == -1, "empty string should have index -1");
    require(empty.pathway.empty(), "empty string should have no pathway");

    const Result singleton = calculate("x");
    require(singleton.assemblyIndex == 0, "one symbol should have index 0");
    require(singleton.pathway.empty(), "one symbol should have no pathway");

    const Result unique = calculate("abcdef");
    require(unique.assemblyIndex == 5, "six unique symbols should have index 5");
    require(unique.pathway.empty(), "unique symbols should have no pathway");

    const Result pair = calculate("abab");
    require(pair.assemblyIndex == 2, "two copies of 'ab' should have index 2");
    requireConsistentPathway("abab", pair);

    require(
        calculate("aaa").assemblyIndex == 2,
        "overlapping copies of 'aa' must not count as a duplication"
    );
    require(
        calculate("aaaa").assemblyIndex == 2,
        "adjacent copies of 'aa' should count as a duplication"
    );
}

void requireMatchesReference(
    const std::string &input,
    ReferenceSearch &orientationSensitiveReference,
    ReferenceSearch &orientationIndependentReference
)
{
    const Result orientationSensitive = calculate(input);
    const int expectedSensitive =
        orientationSensitiveReference.assemblyIndex(input);
    require(
        orientationSensitive.assemblyIndex == expectedSensitive,
        "reference disagreement for '" + input + "': expected " +
            std::to_string(expectedSensitive) + ", got " +
            std::to_string(orientationSensitive.assemblyIndex)
    );
    require(
        !orientationSensitive.runtimeLimitReached &&
            !orientationSensitive.interrupted,
        "unlimited reference case stopped for '" + input + "'"
    );
    requireConsistentPathway(input, orientationSensitive);

    Options options;
    options.acceptReversed = true;
    const Result orientationIndependent = calculate(input, options);
    const int expectedIndependent =
        orientationIndependentReference.assemblyIndex(input);
    require(
        orientationIndependent.assemblyIndex == expectedIndependent,
        "reversal-aware reference disagreement for '" + input +
            "': expected " + std::to_string(expectedIndependent) + ", got " +
            std::to_string(orientationIndependent.assemblyIndex)
    );
    require(
        orientationIndependent.assemblyIndex <=
            orientationSensitive.assemblyIndex,
        "acceptReversed increased the index for '" + input + "'"
    );
    require(
        !orientationIndependent.runtimeLimitReached &&
            !orientationIndependent.interrupted,
        "unlimited reversal-aware case stopped for '" + input + "'"
    );
    requireConsistentPathway(input, orientationIndependent, true);
}

void testExhaustiveShortStrings()
{
    ReferenceSearch orientationSensitiveReference(false);
    ReferenceSearch orientationIndependentReference(true);

    const auto checkAlphabet = [&] (
        const std::string &alphabet,
        size_t maximumLength
    )
    {
        size_t combinationCount = 1;
        for (size_t length = 0; length <= maximumLength; length++)
        {
            if (length > 0) combinationCount *= alphabet.size();
            for (size_t encoded = 0; encoded < combinationCount; encoded++)
            {
                size_t remaining = encoded;
                std::string input(length, alphabet.front());
                for (char &symbol : input)
                {
                    symbol = alphabet[remaining % alphabet.size()];
                    remaining /= alphabet.size();
                }
                requireMatchesReference(
                    input,
                    orientationSensitiveReference,
                    orientationIndependentReference
                );
            }
        }
    };

    checkAlphabet("ab", 8);
    checkAlphabet("abc", 6);
}

void testIntervalUtilities()
{
    implementation::FixedIntervalMap intervals(3);
    require(intervals.empty(), "a new fixed interval map should be empty");
    intervals.insert(5);
    require(!intervals.empty(), "an insertion did not populate its interval map");
    require(
        !intervals.containsTwoDisjointInsertions(),
        "one insertion was mistaken for two disjoint insertions"
    );
    intervals.insert(11);
    require(
        intervals.containsTwoDisjointInsertions(),
        "two separated insertions were not detected"
    );
    requireIntervals(
        intervals.intervals(),
        {{5, 3}, {11, 3}},
        "separated insertion intervals"
    );
    intervals.insert(8);
    intervals.insert(6);
    requireIntervals(
        intervals.intervals(),
        {{5, 9}},
        "touching and overlapping insertion intervals"
    );

    implementation::FixedIntervalMap overlapping(3);
    overlapping.insert(5);
    overlapping.insert(6);
    overlapping.insert(7);
    require(
        !overlapping.containsTwoDisjointInsertions(),
        "overlapping insertions were counted by their multiplicity"
    );
    overlapping.insert(8);
    require(
        overlapping.containsTwoDisjointInsertions(),
        "a six-symbol union did not contain two length-three insertions"
    );

    implementation::FixedIntervalMap singletonRuns(1);
    singletonRuns.insert(7);
    singletonRuns.insert(3);
    singletonRuns.insert(2);
    singletonRuns.insert(4);
    requireIntervals(
        singletonRuns.intervals(),
        {{2, 3}, {7, 1}},
        "out-of-order singleton insertion runs"
    );
    requireIntervals(
        singletonRuns.intervals(true),
        {{2, 3}},
        "singleton run filtering"
    );

    const auto removed = [] (int offset, int length)
    {
        return PathwayStep{{0, length}, {offset, length}};
    };
    requireIntervals(
        implementation::remnantIntervals(10, {}),
        {{0, 10}},
        "remnants without removals"
    );
    requireIntervals(
        implementation::remnantIntervals(
            10,
            {removed(6, 2), removed(2, 2)}
        ),
        {{0, 2}, {4, 2}, {8, 2}},
        "remnants from unsorted separated removals"
    );
    requireIntervals(
        implementation::remnantIntervals(
            10,
            {removed(2, 2), removed(4, 3)}
        ),
        {{0, 2}, {7, 3}},
        "remnants from adjacent removals"
    );
    requireIntervals(
        implementation::remnantIntervals(
            10,
            {removed(4, 4), removed(2, 4), removed(3, 2)}
        ),
        {{0, 2}, {8, 2}},
        "remnants from overlapping and nested removals"
    );
    requireIntervals(
        implementation::remnantIntervals(
            10,
            {removed(8, 2), removed(0, 2)}
        ),
        {{2, 6}},
        "remnants after prefix and suffix removals"
    );
    requireIntervals(
        implementation::remnantIntervals(10, {removed(0, 10)}),
        {},
        "remnants after a full removal"
    );
}

void testMultiStepPathways()
{
    struct TestCase
    {
        std::string input;
        int expectedIndex;
    };
    const std::vector<TestCase> cases{
        {"ababcdcd", 5},
        {"abcababc", 4}
    };

    for (const TestCase &testCase : cases)
    {
        const Result result = calculate(testCase.input);
        require(
            result.assemblyIndex == testCase.expectedIndex,
            "multi-step case '" + testCase.input + "' has index " +
                std::to_string(result.assemblyIndex)
        );
        require(
            result.pathway.size() == 2,
            "multi-step case '" + testCase.input +
                "' should contain two pathway steps"
        );
        requireConsistentPathway(testCase.input, result);
        size_t remnantSymbols = 0;
        for (
            const Interval &interval : implementation::remnantIntervals(
                testCase.input.size(),
                result.pathway
            )
        ) remnantSymbols += static_cast<size_t>(interval.length);
        size_t removedSymbols = 0;
        for (const PathwayStep &step : result.pathway)
            removedSymbols += static_cast<size_t>(step.duplicate.length);
        require(
            remnantSymbols + removedSymbols == testCase.input.size(),
            "multi-step case '" + testCase.input +
                "' lost symbols when generating remnants"
        );
    }
}

void testReverseEquivalence()
{
    const std::string input = "abcabzzxyyxabcab";
    const std::string reversed(input.rbegin(), input.rend());
    const Result forward = calculate(input);
    const Result backward = calculate(reversed);
    require(
        forward.assemblyIndex == backward.assemblyIndex,
        "reversing the whole string changed its assembly index"
    );
    requireConsistentPathway(input, forward);
    requireConsistentPathway(reversed, backward);

    const std::string reverseOnly = "abcxcba";
    const Result orientationSensitive = calculate(reverseOnly);
    Options options;
    options.acceptReversed = true;
    const Result orientationIndependent = calculate(reverseOnly, options);
    require(
        orientationSensitive.assemblyIndex == 6,
        "orientation-sensitive index for 'abcxcba' should be 6"
    );
    require(
        orientationIndependent.assemblyIndex == 4,
        "reversal equivalence should lower the 'abcxcba' index to 4"
    );
    requireConsistentPathway(reverseOnly, orientationIndependent, true);

    const Result reversedAgain = calculate(
        std::string(reverseOnly.rbegin(), reverseOnly.rend()),
        options
    );
    require(
        orientationIndependent.assemblyIndex == reversedAgain.assemblyIndex,
        "acceptReversed result depends on whole-string orientation"
    );
}

bool cancelImmediately()
{
    return true;
}

int cancellationPolls = 0;
int cancellationPollLimit = 0;

bool cancelAfterConfiguredPolls()
{
    cancellationPolls++;
    return cancellationPolls == cancellationPollLimit;
}

void testSearchStops()
{
    const std::string input = numberedBlocks(25);

    Options runtimeOptions;
    runtimeOptions.runtimeTicks = 0;
    const Result timedOut = calculate(input, runtimeOptions);
    require(timedOut.runtimeLimitReached, "zero runtime limit was ignored");
    require(!timedOut.interrupted, "runtime limit reported cancellation");
    require(
        timedOut.assemblyIndex == static_cast<int>(input.size()) - 1,
        "zero-runtime result should retain the initial upper bound"
    );
    require(timedOut.pathway.empty(), "zero-runtime search returned a pathway");
    for (const std::string &noDuplicates : {std::string(), std::string("abcdef")})
    {
        const Result stopped = calculate(noDuplicates, runtimeOptions);
        require(
            stopped.runtimeLimitReached,
            "zero runtime was ignored without duplicate substrings"
        );
        require(!stopped.interrupted, "zero runtime reported cancellation");
        require(
            stopped.assemblyIndex ==
                static_cast<int>(noDuplicates.size()) - 1,
            "zero-runtime no-duplicate result changed its initial bound"
        );
        require(
            stopped.pathway.empty(),
            "zero-runtime no-duplicate result returned a pathway"
        );
    }

    Options cancellationOptions;
    cancellationOptions.cancellationRequested = &cancelImmediately;
    const Result cancelled = calculate(input, cancellationOptions);
    require(cancelled.interrupted, "cancellation callback was ignored");
    require(!cancelled.runtimeLimitReached, "cancellation reported a timeout");
    require(
        cancelled.assemblyIndex == static_cast<int>(input.size()) - 1,
        "immediately cancelled result should retain the initial upper bound"
    );
    require(
        cancelled.pathway.empty(),
        "immediately cancelled search returned a pathway"
    );
    for (const std::string &noDuplicates : {std::string(), std::string("abcdef")})
    {
        const Result stopped = calculate(noDuplicates, cancellationOptions);
        require(
            stopped.interrupted,
            "cancellation was ignored without duplicate substrings"
        );
        require(
            !stopped.runtimeLimitReached,
            "no-duplicate cancellation reported a timeout"
        );
        require(
            stopped.assemblyIndex ==
                static_cast<int>(noDuplicates.size()) - 1,
            "cancelled no-duplicate result changed its initial bound"
        );
        require(
            stopped.pathway.empty(),
            "cancelled no-duplicate result returned a pathway"
        );
    }

    Options competingOptions;
    competingOptions.runtimeTicks = 0;
    competingOptions.cancellationRequested = &cancelImmediately;
    const Result cancellationWins = calculate(input, competingOptions);
    require(
        cancellationWins.interrupted,
        "cancellation was not reported when two stop conditions applied"
    );
    require(
        !cancellationWins.runtimeLimitReached,
        "a lower-priority runtime limit obscured cancellation"
    );

    cancellationPolls = 0;
    cancellationPollLimit = 2;
    Options postPreprocessCancellationOptions;
    postPreprocessCancellationOptions.cancellationRequested =
        &cancelAfterConfiguredPolls;
    const Result stoppedAfterPreprocess = calculate(
        "abcdef",
        postPreprocessCancellationOptions
    );
    require(
        cancellationPolls == cancellationPollLimit &&
            stoppedAfterPreprocess.interrupted,
        "a no-duplicate input was not polled after preprocessing"
    );
    require(
        stoppedAfterPreprocess.assemblyIndex == 5 &&
            stoppedAfterPreprocess.pathway.empty(),
        "post-preprocessing cancellation changed a no-duplicate result"
    );

    cancellationPolls = 0;
    cancellationPollLimit = 4;
    Options delayedCancellationOptions;
    delayedCancellationOptions.cancellationRequested =
        &cancelAfterConfiguredPolls;
    const std::string partialInput = "abcabcabc";
    const Result partiallySearched = calculate(
        partialInput,
        delayedCancellationOptions
    );
    require(
        cancellationPolls == cancellationPollLimit,
        "delayed cancellation used an unexpected number of polls"
    );
    require(partiallySearched.interrupted, "delayed cancellation was ignored");
    require(
        !partiallySearched.runtimeLimitReached,
        "delayed cancellation reported a timeout"
    );
    require(
        partiallySearched.assemblyIndex >= 4 &&
            partiallySearched.assemblyIndex <
                static_cast<int>(partialInput.size()),
        "delayed cancellation returned an invalid best-so-far index"
    );
    requireConsistentPathway(partialInput, partiallySearched);
}

void testJsonAndPathwayOutput()
{
    const std::vector<std::string> controlEscapes{
        "\\u0000", "\\u0001", "\\u0002", "\\u0003",
        "\\u0004", "\\u0005", "\\u0006", "\\u0007",
        "\\b", "\\t", "\\n", "\\u000B", "\\f", "\\r",
        "\\u000E", "\\u000F", "\\u0010", "\\u0011",
        "\\u0012", "\\u0013", "\\u0014", "\\u0015",
        "\\u0016", "\\u0017", "\\u0018", "\\u0019",
        "\\u001A", "\\u001B", "\\u001C", "\\u001D",
        "\\u001E", "\\u001F"
    };
    for (size_t value = 0; value < controlEscapes.size(); value++)
    {
        std::ostringstream controlOutput;
        implementation::writeJsonString(
            std::string(1, static_cast<char>(value)),
            controlOutput
        );
        require(
            controlOutput.str() == '"' + controlEscapes[value] + '"',
            "incorrect JSON escape for control byte " +
                std::to_string(value) + ": " + controlOutput.str()
        );
    }

    std::string special = "quote\"\\\b\f\n\r\t";
    special.push_back('\x01');
    special.push_back('\x1f');
    special += " end";

    std::ostringstream encoded;
    implementation::writeJsonString(special, encoded);
    require(
        encoded.str() ==
            "\"quote\\\"\\\\\\b\\f\\n\\r\\t\\u0001\\u001F end\"",
        "JSON string escaping is incorrect: " + encoded.str()
    );

    const std::string input = special + '|' + special;
    const Result result = calculate(input);
    requireConsistentPathway(input, result);

    const auto nonce = std::chrono::steady_clock::now()
        .time_since_epoch().count();
    const std::filesystem::path outputPath =
        std::filesystem::temp_directory_path() /
        ("assemblycpp-string-pathway-" + std::to_string(nonce) + ".json");
    const std::filesystem::path multiStepOutputPath =
        std::filesystem::temp_directory_path() /
        ("assemblycpp-string-pathway-multi-step-" +
            std::to_string(nonce) + ".json");

    struct RemoveFile
    {
        std::filesystem::path path;
        ~RemoveFile()
        {
            std::error_code ignored;
            std::filesystem::remove(path, ignored);
        }
    } cleanup{outputPath}, multiStepCleanup{multiStepOutputPath};

    std::string error;
    require(
        writePathway(outputPath.string(), input, result, error),
        "could not write pathway JSON: " + error
    );

    std::ifstream pathwayFile(outputPath);
    require(pathwayFile.is_open(), "pathway JSON file was not created");
    std::ostringstream contents;
    contents << pathwayFile.rdbuf();
    require(pathwayFile.good() || pathwayFile.eof(), "could not read pathway JSON");

    std::ostringstream encodedInput;
    implementation::writeJsonString(input, encodedInput);
    const std::string json = contents.str();
    require(
        json.find(encodedInput.str()) != std::string::npos,
        "pathway JSON does not contain the correctly escaped input"
    );
    require(
        json.find("\"remnant\"") != std::string::npos &&
            json.find("\"duplicates\"") != std::string::npos,
        "pathway JSON is missing required sections"
    );

    Result synthetic;
    synthetic.assemblyIndex = 6;
    synthetic.pathway = {
        {{0, 2}, {2, 2}},
        {{4, 3}, {7, 3}}
    };
    requireConsistentPathway("ababXYZXYZ", synthetic);
    error.clear();
    require(
        writePathway(
            multiStepOutputPath.string(),
            "ababXYZXYZ",
            synthetic,
            error
        ),
        "could not write multi-step pathway fixture: " + error
    );
    std::ifstream multiStepFile(multiStepOutputPath);
    require(
        multiStepFile.is_open(),
        "multi-step pathway JSON file was not created"
    );
    std::ostringstream multiStepContents;
    multiStepContents << multiStepFile.rdbuf();
    require(
        multiStepFile.good() || multiStepFile.eof(),
        "could not read multi-step pathway JSON"
    );
    const std::string compact = withoutJsonFormatting(multiStepContents.str());
    require(
        compact.find("\"Fragments\":[\"ababXYZXYZ\"]") !=
            std::string::npos,
        "multi-step pathway JSON changed its input fragment"
    );
    require(
        compact.find("\"Fragments\":[\"ab\",\"XYZ\"]") !=
            std::string::npos &&
            compact.find("\"Positions\":[0,4]") != std::string::npos,
        "multi-step pathway JSON changed its remnant fragments"
    );
    size_t firstDuplicate = compact.find(
        "{\"Left\":[0,2],\"Right\":[2,2]}"
    );
    if (firstDuplicate == std::string::npos)
        firstDuplicate = compact.find(
            "{\"Right\":[2,2],\"Left\":[0,2]}"
        );
    size_t secondDuplicate = compact.find(
        "{\"Left\":[4,3],\"Right\":[7,3]}"
    );
    if (secondDuplicate == std::string::npos)
        secondDuplicate = compact.find(
            "{\"Right\":[7,3],\"Left\":[4,3]}"
        );
    require(
        firstDuplicate != std::string::npos &&
            secondDuplicate != std::string::npos &&
            firstDuplicate < secondDuplicate,
        "multi-step pathway JSON changed its ordered duplicates"
    );

    const std::filesystem::path missingDirectory =
        std::filesystem::temp_directory_path() /
        ("assemblycpp-string-pathway-missing-" + std::to_string(nonce));
    require(
        !std::filesystem::exists(missingDirectory),
        "temporary missing-directory fixture unexpectedly exists"
    );
    const std::filesystem::path unavailablePath =
        missingDirectory / "pathway.json";
    std::string openError;
    require(
        !writePathway(
            unavailablePath.string(),
            "abab",
            calculate("abab"),
            openError
        ),
        "pathway writing unexpectedly created a missing parent directory"
    );
    require(
        openError ==
            "could not open output file '" + unavailablePath.string() + "'",
        "pathway open failure returned the wrong error: " + openError
    );
}

} // namespace

int main()
{
    try
    {
        testUpstreamRegressionCases();
        testTrivialStrings();
        testExhaustiveShortStrings();
        testIntervalUtilities();
        testMultiStepPathways();
        testReverseEquivalence();
        testSearchStops();
        testJsonAndPathwayOutput();
    }
    catch (const std::exception &error)
    {
        std::cerr << "string assembly test failed: " << error.what() << '\n';
        return 1;
    }

    std::cout << "string assembly tests passed\n";
    return 0;
}
