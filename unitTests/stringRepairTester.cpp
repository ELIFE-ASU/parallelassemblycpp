#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../src/stringRepair.h"

namespace repair = parallelassemblycpp::detail::stringRepair;

namespace
{

void require(bool condition, const std::string &message)
{
    if (!condition) throw std::runtime_error(message);
}

/** Independent exhaustive duplicate-removal oracle for small ASCII inputs.
 * Every state is a multiset of literal fragments. Removing two equal substrings
 * retains one copy and saves length-1 joins; all non-overlapping pairs are tried.
 * This deliberately uses neither production string search nor Re-Pair rules.
 */
class Oracle
{
    struct Occurrence { size_t fragment, start; };
    bool reversed_;
    std::unordered_map<std::string, int> memo_;

    std::string canonical(std::string text) const
    {
        if (reversed_)
            text = std::min(text, std::string(text.rbegin(), text.rend()));
        return text;
    }

    std::string key(const std::vector<std::string> &fragments) const
    {
        std::vector<std::string> sorted;
        for (const auto &fragment : fragments)
            if (fragment.size() > 1) sorted.push_back(canonical(fragment));
        std::sort(sorted.begin(), sorted.end());
        std::string result;
        for (const auto &fragment : sorted)
            result += std::to_string(fragment.size()) + ':' + fragment + ';';
        return result;
    }

    int savings(const std::vector<std::string> &fragments)
    {
        const std::string state = key(fragments);
        const auto known = memo_.find(state);
        if (known != memo_.end()) return known->second;
        int best = 0;
        size_t maximum = 0;
        for (const auto &fragment : fragments)
            maximum = std::max(maximum, fragment.size());
        std::unordered_set<std::string> visited;
        for (size_t length = 2; length <= maximum; ++length)
        {
            std::vector<Occurrence> occurrences;
            for (size_t f = 0; f < fragments.size(); ++f)
                for (size_t start = 0; start + length <= fragments[f].size(); ++start)
                    occurrences.push_back({f, start});
            for (size_t a = 0; a < occurrences.size(); ++a)
                for (size_t b = a + 1; b < occurrences.size(); ++b)
                {
                    const auto left = occurrences[a];
                    const auto right = occurrences[b];
                    if (left.fragment == right.fragment && left.start + length > right.start)
                        continue;
                    const auto text = fragments[left.fragment].substr(left.start, length);
                    if (canonical(text) != canonical(
                            fragments[right.fragment].substr(right.start, length)))
                        continue;
                    std::vector<std::string> next{text};
                    for (size_t f = 0; f < fragments.size(); ++f)
                    {
                        std::vector<size_t> cuts;
                        if (f == left.fragment) cuts.push_back(left.start);
                        if (f == right.fragment) cuts.push_back(right.start);
                        size_t start = 0;
                        for (const size_t cut : cuts)
                        {
                            if (cut - start > 1)
                                next.push_back(fragments[f].substr(start, cut - start));
                            start = cut + length;
                        }
                        if (fragments[f].size() - start > 1)
                            next.push_back(fragments[f].substr(start));
                    }
                    const auto transition = std::to_string(length) + ':' + key(next);
                    if (!visited.insert(transition).second) continue;
                    best = std::max(best, static_cast<int>(length) - 1 + savings(next));
                }
        }
        memo_.emplace(state, best);
        return best;
    }

public:
    explicit Oracle(bool acceptReversed) : reversed_(acceptReversed) {}
    int index(const std::string &text)
    {
        return static_cast<int>(text.size()) - 1 - savings({text});
    }
};

std::u32string oriented(std::u32string text, bool reversed)
{
    if (reversed) std::reverse(text.begin(), text.end());
    return text;
}

std::string encode(const std::u32string &text)
{
    std::string bytes;
    for (const char32_t cp : text)
    {
        if (cp <= 0x7f) bytes.push_back(static_cast<char>(cp));
        else if (cp <= 0x7ff)
        {
            bytes.push_back(static_cast<char>(0xc0 | (cp >> 6)));
            bytes.push_back(static_cast<char>(0x80 | (cp & 0x3f)));
        }
        else if (cp <= 0xffff)
        {
            bytes.push_back(static_cast<char>(0xe0 | (cp >> 12)));
            bytes.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3f)));
            bytes.push_back(static_cast<char>(0x80 | (cp & 0x3f)));
        }
        else
        {
            bytes.push_back(static_cast<char>(0xf0 | (cp >> 18)));
            bytes.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3f)));
            bytes.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3f)));
            bytes.push_back(static_cast<char>(0x80 | (cp & 0x3f)));
        }
    }
    return bytes;
}

void replay(const std::u32string &input, const repair::Result &result, bool reversed)
{
    require(result.acceptReversed == reversed, "certificate changed reversal mode");
    require(result.trivialUpperBound == static_cast<int>(input.size()) - 1,
        "trivial bound counts bytes rather than Unicode scalars");
    require(result.ruleCount == static_cast<int>(result.rules.size()), "wrong rule count");
    require(result.remainingFragments == static_cast<int>(result.residual.size()),
        "wrong residual count");
    require(result.upperBound == result.ruleCount + result.remainingFragments - 1,
        "bound disagrees with grammar construction cost");
    require(result.upperBound <= result.trivialUpperBound, "bound exceeds trivial construction");
    std::map<int, std::u32string> symbols;
    std::u32string firstAppearance;
    for (const char32_t cp : input)
        if (firstAppearance.find(cp) == std::u32string::npos) firstAppearance.push_back(cp);
    require(firstAppearance.size() == result.terminals.size(), "wrong terminal alphabet");
    for (size_t i = 0; i < result.terminals.size(); ++i)
    {
        const auto &terminal = result.terminals[i];
        require(terminal.codePoint == firstAppearance[i], "terminal order is not first appearance");
        require(symbols.emplace(terminal.id, std::u32string(1, terminal.codePoint)).second,
            "duplicate terminal id");
    }
    for (const auto &rule : result.rules)
    {
        require(symbols.count(rule.left) && symbols.count(rule.right),
            "rule refers to a missing or forward symbol");
        require(reversed || (!rule.leftReversed && !rule.rightReversed),
            "directed rule uses a reversed child");
        const auto expansion = oriented(symbols.at(rule.left), rule.leftReversed) +
            oriented(symbols.at(rule.right), rule.rightReversed);
        require(static_cast<int>(expansion.size()) == rule.length, "incorrect rule length");
        require(rule.length > 1, "rule does not join two nonempty fragments");
        require(!reversed || expansion <= oriented(expansion, true),
            "reverse-equivalent rule has a noncanonical expansion");
        require(symbols.emplace(rule.id, expansion).second, "duplicate rule id");
    }
    std::u32string reconstructed;
    for (const auto &fragment : result.residual)
    {
        require(symbols.count(fragment.symbol), "residual refers to missing symbol");
        require(reversed || !fragment.reversed, "directed residual reverses a symbol");
        const auto expansion = oriented(symbols.at(fragment.symbol), fragment.reversed);
        require(fragment.offset == static_cast<int>(reconstructed.size()),
            "residual intervals are not contiguous scalar offsets");
        require(fragment.length == static_cast<int>(expansion.size()),
            "residual length differs from grammar expansion");
        require(input.substr(reconstructed.size(), expansion.size()) == expansion,
            "residual expansion differs from original substring");
        reconstructed += expansion;
    }
    require(reconstructed == input, "certificate does not reconstruct the original input");
}

std::string json(const repair::Result &result, const std::string &input)
{
    std::ostringstream output;
    repair::writeJson(result, input, output);
    require(output.good(), "JSON writer failed");
    return output.str();
}

void testKnownCases()
{
    const std::vector<std::pair<std::string, int>> cases{
        {"", -1}, {"a", 0}, {"ab", 1}, {"aa", 1}, {"aaa", 2},
        {"aaaa", 2}, {"aaaaa", 3}, {"abab", 2}, {"ababa", 3},
        {"abcabc", 3}, {"aaaaaaaa", 3}, {"abcdef", 5},
        {"abababababababab", 4}, {std::string(16, 'a'), 4},
        {std::string(17, 'a'), 5}
    };
    for (const auto &[input, expected] : cases)
        for (const bool reversed : {false, true})
        {
            const auto result = repair::calculate(input, reversed);
            require(result.upperBound == expected, "wrong bound for known input '" + input + "'");
            replay(std::u32string(input.begin(), input.end()), result, reversed);
        }
    require(repair::calculate("abba").upperBound == 3, "directed search reused a reverse copy");
    require(repair::calculate("abba", true).upperBound == 2, "reversed pair was not reused");
    require(repair::calculate("abcxcba", true).upperBound == 4,
        "nested reverse-equivalent rules were not reused");
}

void testExhaustiveBounds()
{
    for (const bool reversed : {false, true})
    {
        Oracle oracle(reversed);
        for (size_t length = 0; length <= 8; ++length)
            for (size_t bits = 0; bits < (size_t{1} << length); ++bits)
            {
                std::string input(length, 'a');
                for (size_t i = 0; i < length; ++i)
                    if (bits & (size_t{1} << i)) input[i] = 'b';
                const auto result = repair::calculate(input, reversed);
                require(result.upperBound >= oracle.index(input),
                    "heuristic returned an impossible bound for '" + input + "'");
                replay(std::u32string(input.begin(), input.end()), result, reversed);
            }
        for (const std::string input : {"abcxcba", "abcabcab", "abcababc", "abccbacb", "ababcdcd"})
        {
            const auto result = repair::calculate(input, reversed);
            require(result.upperBound >= oracle.index(input), "ternary oracle bound violation");
            replay(std::u32string(input.begin(), input.end()), result, reversed);
        }
    }
}

void testUnicodeAndJson()
{
    const std::u32string boundaries{0, 0x7f, 0x80, 0x7ff, 0x800, 0xd7ff,
        0xe000, 0xffff, 0x10000, 0x10ffff};
    const std::vector<std::u32string> cases{
        boundaries, boundaries + boundaries, U"\u00e9a\u00e9a", U"\U0001f600\u00e9\u4e2d\U0001f600\u00e9\u4e2d", U"\u03b1\u03b2\u03b3x\u03b3\u03b2\u03b1",
        U"e\u0301\u00e9", std::u32string{0, 0x1f600, 0, 0x1f600},
        U"a\"\\\b\f\n\r\t\x01z", U"\u03b1\u03b2\u03b1\u03b3\u03b1\u03b3\u03b1\u03b2\u03b1", U"\u03b1\u03b2\u03b2\u03b1"
    };
    for (const auto &input : cases)
        for (const bool reversed : {false, true})
        {
            const auto bytes = encode(input);
            const auto result = repair::calculate(bytes, reversed);
            replay(input, result, reversed);
            const auto certificate = json(result, bytes);
            require(certificate.find("string-repair-assembly-v1") != std::string::npos,
                "JSON omitted the versioned schema");
            require(certificate == json(repair::calculate(bytes, reversed), bytes),
                "JSON certificate is not deterministic");
        }
    require(repair::calculate(encode(boundaries)).upperBound == 9,
        "distinct Unicode scalars were merged");
    require(repair::calculate("e\xcc\x81").upperBound == 1,
        "Unicode combining characters were normalized unexpectedly");
}

void testDeterminismAndLargeInputs()
{
    for (const std::string input : {"ababcdcd", "abccbaabcddcba", "abcababc", "babaabaa"})
        for (const bool reversed : {false, true})
        {
            const auto expected = json(repair::calculate(input, reversed), input);
            for (int iteration = 0; iteration < 12; ++iteration)
                require(json(repair::calculate(input, reversed), input) == expected,
                    "ties produce nondeterministic grammar certificates");
        }
    for (const bool reversed : {false, true})
        for (const bool alternating : {false, true})
        {
            constexpr size_t length = size_t{1} << 18;
            std::string input(length, 'a');
            if (alternating)
                for (size_t i = 1; i < length; i += 2) input[i] = 'b';
            const auto result = repair::calculate(input, reversed);
            require(result.upperBound == 18, "large repeated string did not compress by doubling");
            replay(std::u32string(input.begin(), input.end()), result, reversed);
        }
}

struct NaiveResult
{
    std::vector<std::u32string> rules;
    std::vector<std::u32string> residual;
};

/** Full scans and literal expanded strings make a small, intentionally slow
 * Re-Pair reference. It shares no occurrence/run/heap bookkeeping with the
 * implementation and catches missed profitable replacements as well as ties.
 */
NaiveResult naiveRepair(const std::u32string &input, bool reversed)
{
    auto canonical = [reversed](const std::u32string &text)
    {
        return reversed ? std::min(text, oriented(text, true)) : text;
    };
    NaiveResult result;
    std::map<std::u32string, bool> known;
    for (const char32_t cp : input)
    {
        result.residual.emplace_back(1, cp);
        known[result.residual.back()] = true;
    }
    for (;;)
    {
        struct Candidate
        {
            int count = 0;
            size_t first = 0;
            size_t end = 0;
        };
        std::map<std::u32string, Candidate> candidates;
        for (size_t i = 0; i + 1 < result.residual.size(); ++i)
        {
            const auto text = canonical(result.residual[i] + result.residual[i + 1]);
            auto &candidate = candidates[text];
            if (candidate.count == 0) candidate.first = i;
            if (candidate.end > i) continue;
            ++candidate.count;
            candidate.end = i + 2;
        }
        int bestSaving = 0;
        size_t first = result.residual.size();
        std::u32string winner;
        for (const auto &[text, candidate] : candidates)
        {
            const int saving = candidate.count - (known.count(text) ? 0 : 1);
            if (saving > bestSaving ||
                (saving == bestSaving && saving > 0 && candidate.first < first))
            {
                bestSaving = saving;
                first = candidate.first;
                winner = text;
            }
        }
        if (bestSaving <= 0) break;
        if (known.emplace(winner, true).second) result.rules.push_back(winner);
        std::vector<std::u32string> remaining;
        for (size_t i = 0; i < result.residual.size();)
        {
            if (i + 1 < result.residual.size() &&
                canonical(result.residual[i] + result.residual[i + 1]) == winner)
            {
                remaining.push_back(result.residual[i] + result.residual[i + 1]);
                i += 2;
            }
            else remaining.push_back(result.residual[i++]);
        }
        result.residual = std::move(remaining);
    }
    return result;
}

void compareNaive(const std::u32string &input, bool reversed)
{
    const auto expected = naiveRepair(input, reversed);
    const auto actual = repair::calculate(encode(input), reversed);
    replay(input, actual, reversed);
    require(actual.ruleCount == static_cast<int>(expected.rules.size()),
        "greedy rule count differs from full-scan reference");
    require(actual.remainingFragments == static_cast<int>(expected.residual.size()),
        "greedy residual count differs from full-scan reference");
    std::map<int, std::u32string> expansions;
    for (const auto &terminal : actual.terminals)
        expansions.emplace(terminal.id, std::u32string(1, terminal.codePoint));
    for (size_t i = 0; i < actual.rules.size(); ++i)
    {
        const auto &rule = actual.rules[i];
        const auto text = oriented(expansions.at(rule.left), rule.leftReversed) +
            oriented(expansions.at(rule.right), rule.rightReversed);
        require(text == expected.rules[i], "greedy rule trace differs from full-scan reference");
        expansions.emplace(rule.id, text);
    }
    for (size_t i = 0; i < actual.residual.size(); ++i)
    {
        const auto &fragment = actual.residual[i];
        require(oriented(expansions.at(fragment.symbol), fragment.reversed) == expected.residual[i],
            "greedy residual trace differs from full-scan reference");
    }
}

void testNaiveGreedyReference()
{
    std::uint32_t randomState = 0x914da21U;
    auto random = [&randomState]()
    {
        randomState ^= randomState << 13;
        randomState ^= randomState >> 17;
        randomState ^= randomState << 5;
        return randomState;
    };
    for (const bool reversed : {false, true})
    {
        // The first fixture groups one expansion across unequal child splits;
        // the second profitably reuses an existing rule in a later iteration.
        compareNaive(U"aabbaa", reversed);
        compareNaive(U"cbcacaaccbcacbacbccbccbac", reversed);
        for (size_t length = 0; length <= 8; ++length)
            for (size_t bits = 0; bits < (size_t{1} << length); ++bits)
            {
                std::u32string input(length, U'a');
                for (size_t i = 0; i < length; ++i)
                    if (bits & (size_t{1} << i)) input[i] = U'b';
                compareNaive(input, reversed);
            }
        for (size_t trial = 0; trial < 1000; ++trial)
        {
            const size_t length = random() % 101;
            const char32_t alphabetSize = 1 + random() % 6;
            std::u32string input;
            for (size_t i = 0; i < length; ++i)
                input.push_back(U'a' + random() % alphabetSize);
            compareNaive(input, reversed);
            if (trial % 10 == 0)
                compareNaive(input + oriented(input, true) + input, reversed);
        }
    }
}

void testInvalidUtf8()
{
    const std::vector<std::string> malformed{
        "\x80", "\xbf", "\xc0\xaf", "\xc1\xbf", "\xc2", "\xc2x",
        "\xe0\x80\x80", "\xe0\x9f\xbf", "\xe1\x80", "\xed\xa0\x80",
        "\xed\xbf\xbf", "\xf0\x80\x80\x80", "\xf0\x8f\xbf\xbf",
        "\xf0\x90\x80", "\xf4\x90\x80\x80", "\xf5\x80\x80\x80",
        "\xf8\x88\x80\x80\x80", "\xfe", "\xff", "a\xe2(\xa1z"
    };
    for (const auto &input : malformed)
        for (const bool reversed : {false, true})
        {
            bool rejected = false;
            try { static_cast<void>(repair::calculate(input, reversed)); }
            catch (const std::invalid_argument &) { rejected = true; }
            require(rejected, "malformed UTF-8 was accepted");
        }
}

} // namespace

int main()
{
    try
    {
        testKnownCases();
        testExhaustiveBounds();
        testUnicodeAndJson();
        testDeterminismAndLargeInputs();
        testInvalidUtf8();
        testNaiveGreedyReference();
    }
    catch (const std::exception &error)
    {
        std::cerr << "string Re-Pair test failed: " << error.what() << '\n';
        return 1;
    }
    std::cout << "string Re-Pair tests passed\n";
}
