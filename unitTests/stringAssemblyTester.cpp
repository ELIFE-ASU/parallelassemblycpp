#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../v5/stringAssembly.h"

using assemblycpp::detail::stringAssembly::Interval;
using assemblycpp::detail::stringAssembly::Options;
using assemblycpp::detail::stringAssembly::PathwayStep;
using assemblycpp::detail::stringAssembly::Result;
using assemblycpp::detail::stringAssembly::calculate;
using assemblycpp::detail::stringAssembly::writePathway;

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

        const std::string match = intervalText(input, step.match);
        const std::string duplicate = intervalText(input, step.duplicate);
        std::string reversedDuplicate(duplicate.rbegin(), duplicate.rend());
        require(
            match == duplicate || (acceptReversed && match == reversedDuplicate),
            prefix + "does not identify equivalent text"
        );
        duplicatedSymbols += step.match.length - 1;
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

    Options cancellationOptions;
    cancellationOptions.cancellationRequested = &cancelImmediately;
    const Result cancelled = calculate(input, cancellationOptions);
    require(cancelled.interrupted, "cancellation callback was ignored");
    require(!cancelled.runtimeLimitReached, "cancellation reported a timeout");
    require(
        cancelled.assemblyIndex == static_cast<int>(input.size()) - 1,
        "immediately cancelled result should retain the initial upper bound"
    );
    require(cancelled.pathway.empty(), "immediately cancelled search returned a pathway");
}

void testJsonAndPathwayOutput()
{
    std::string special = "quote\"\\\b\f\n\r\t";
    special.push_back('\x01');
    special.push_back('\x1f');
    special += " end";

    std::ostringstream encoded;
    assemblycpp::detail::stringAssembly::implementation::writeJsonString(
        special,
        encoded
    );
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

    struct RemoveFile
    {
        std::filesystem::path path;
        ~RemoveFile()
        {
            std::error_code ignored;
            std::filesystem::remove(path, ignored);
        }
    } cleanup{outputPath};

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
    assemblycpp::detail::stringAssembly::implementation::writeJsonString(
        input,
        encodedInput
    );
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
}

} // namespace

int main()
{
    try
    {
        testUpstreamRegressionCases();
        testTrivialStrings();
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
