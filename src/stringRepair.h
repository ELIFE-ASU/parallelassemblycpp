#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "stringEncoding.h"

/**
 * Constructive Re-Pair assembly bounds for ordered Unicode scalar strings.
 *
 * Adjacent active fragments are grouped by their FULL expanded text, so the
 * same product can be reused even when it has different binary decompositions.
 * With acceptReversed, a fragment and its scalar reversal share one symbol.
 * Each new binary production costs one join; an existing product costs zero.
 * We choose the product giving the largest immediate saving, breaking ties by
 * its first occurrence in the input, and replace a maximum, leftmost set of
 * nonoverlapping occurrences. The remaining ordered fragments need r-1 joins.
 * Consequently rules+r-1 is a witnessed upper bound, not an exact minimum.
 * The empty string has bound -1, matching the exact string solver.
 *
 * Collision-free doubling ranks identify arbitrary substrings in constant
 * time after O(n log n) preprocessing. Occurrence runs, an active-position
 * Fenwick tree, and a priority queue maintain exact nonoverlap frequencies
 * using local updates. Total time and space are O(n log n), including the
 * rank table; this is not a claim of the linear-space Re-Pair algorithm.
 * All state belongs to one calculation; independent calls are thread-safe.
 */
namespace parallelassemblycpp::detail::stringRepair
{
struct Terminal
{
    int id = -1;
    char32_t codePoint = 0;
};

/** One join, with each child optionally reversed before concatenation. */
struct Rule
{
    int id = -1;
    int left = -1;
    int right = -1;
    bool leftReversed = false;
    bool rightReversed = false;
    int length = 0;
};

/** A half-open input interval in Unicode scalar positions, not UTF-8 bytes. */
struct Fragment
{
    int symbol = -1;
    bool reversed = false;
    int offset = 0;
    int length = 0;
};

struct Result
{
    int upperBound = -1;
    int trivialUpperBound = -1;
    int ruleCount = 0;
    int remainingFragments = 0;
    bool acceptReversed = false;
    std::vector<Terminal> terminals;
    std::vector<Rule> rules;
    std::vector<Fragment> residual;
};

namespace implementation
{
/** An interrupted construction has no complete certificate to publish. */
struct Cancelled : std::exception
{
    const char *what() const noexcept override
    {
        return "string Re-Pair construction cancelled";
    }
};

/** Amortize callback costs while bounding work between cooperative polls. */
class Cancellation
{
    std::function<bool()> requested_;
    mutable std::size_t operations_ = 0;

public:
    explicit Cancellation(std::function<bool()> requested):
        requested_(std::move(requested)) {}

    void poll() const
    {
        if (requested_ && requested_()) throw Cancelled{};
    }

    void tick() const
    {
        if (requested_ && (++operations_ & 1023U) == 0) poll();
    }
};

/** Two overlapping power-of-two blocks exactly cover a substring. */
struct Key
{
    int length = 0;
    std::uint32_t first = 0;
    std::uint32_t last = 0;

    bool operator==(const Key &) const = default;
    bool operator<(const Key &other) const
    {
        return std::tie(length, first, last) <
            std::tie(other.length, other.first, other.last);
    }
};

class Substrings
{
    std::size_t length_;
    bool acceptReversed_;
    std::vector<unsigned char> logarithms_;
    std::vector<std::vector<std::uint32_t>> ranks_;

    Key forward(std::size_t offset, int length) const
    {
        const auto level = logarithms_[static_cast<std::size_t>(length)];
        const std::size_t block = std::size_t{1} << level;
        const auto &rank = ranks_[level];
        return {length, rank[offset],
                rank[offset + static_cast<std::size_t>(length) - block]};
    }

public:
    Substrings(const std::u32string &input, bool acceptReversed,
               const Cancellation &cancellation):
        length_(input.size()), acceptReversed_(acceptReversed),
        logarithms_(length_ + 1, 0)
    {
        cancellation.poll();
        for (std::size_t i = 2; i <= length_; ++i)
        {
            cancellation.tick();
            logarithms_[i] = static_cast<unsigned char>(logarithms_[i / 2] + 1);
        }
        if (input.empty()) return;

        // The two halves need no separator: every queried interval remains
        // inside its own half. uint32 ranks suffice for 2*(INT_MAX-1) scalars.
        const std::size_t count = length_ * (acceptReversed ? 2U : 1U);
        const auto scalar = [&](std::size_t position) {
            return position < length_ ? input[position] :
                input[count - position - 1];
        };
        std::vector<std::size_t> order(count), temporary(count);
        for (std::size_t i = 0; i < count; ++i)
        {
            cancellation.tick();
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            cancellation.tick();
            return scalar(a) < scalar(b);
        });
        cancellation.poll();
        ranks_.emplace_back(count);
        auto &initial = ranks_.back();
        std::uint32_t classes = 0;
        char32_t previous = 0;
        for (std::size_t i = 0; i < count; ++i)
        {
            cancellation.tick();
            const char32_t value = scalar(order[i]);
            if (i == 0 || value != previous) ++classes;
            initial[order[i]] = classes;
            previous = value;
        }

        std::vector<std::size_t> buckets;
        for (std::size_t block = 1; block <= length_ / 2; block *= 2)
        {
            cancellation.poll();
            const auto &old = ranks_.back();
            const auto second = [&](std::size_t i) -> std::uint32_t {
                return i + block < count ? old[i + block] : 0;
            };
            // Stable counting sorts of the second and then first rank keep
            // each doubling level linear, without probabilistic hashing.
            buckets.assign(static_cast<std::size_t>(classes) + 1, 0);
            for (std::size_t i = 0; i < count; ++i)
            {
                cancellation.tick();
                ++buckets[second(i)];
            }
            std::size_t total = 0;
            for (auto &bucket : buckets)
            {
                cancellation.tick();
                const auto frequency = bucket;
                bucket = total;
                total += frequency;
            }
            for (std::size_t i = 0; i < count; ++i)
            {
                cancellation.tick();
                temporary[buckets[second(i)]++] = i;
            }
            for (auto &bucket : buckets)
            {
                cancellation.tick();
                bucket = 0;
            }
            for (const auto i : temporary)
            {
                cancellation.tick();
                ++buckets[old[i]];
            }
            total = 0;
            for (auto &bucket : buckets)
            {
                cancellation.tick();
                const auto frequency = bucket;
                bucket = total;
                total += frequency;
            }
            for (const auto i : temporary)
            {
                cancellation.tick();
                order[buckets[old[i]]++] = i;
            }

            std::vector<std::uint32_t> next(count);
            classes = 0;
            for (std::size_t i = 0; i < count; ++i)
            {
                cancellation.tick();
                if (i == 0 || old[order[i]] != old[order[i - 1]] ||
                    second(order[i]) != second(order[i - 1])) ++classes;
                next[order[i]] = classes;
            }
            ranks_.push_back(std::move(next));
        }
    }

    std::pair<Key, bool> canonical(int offset, int length) const
    {
        const Key direct = forward(static_cast<std::size_t>(offset), length);
        if (!acceptReversed_) return {direct, false};
        const Key reverse = forward(2 * length_ -
            static_cast<std::size_t>(offset) - static_cast<std::size_t>(length), length);
        return reverse < direct ? std::pair{reverse, true} :
            std::pair{direct, false};
    }

    bool palindrome(const Fragment &fragment) const
    {
        if (fragment.length == 1) return true;
        if (!acceptReversed_) return false;
        return forward(static_cast<std::size_t>(fragment.offset), fragment.length) ==
            forward(2 * length_ - static_cast<std::size_t>(fragment.offset) -
                static_cast<std::size_t>(fragment.length), fragment.length);
    }
};

/** Number of live token starts in any original-position interval. */
class ActivePositions
{
    std::vector<int> tree_;

    int prefix(int end) const
    {
        int count = 0;
        for (std::size_t i = static_cast<std::size_t>(end); i != 0; i &= i - 1)
            count += tree_[i];
        return count;
    }

public:
    ActivePositions(std::size_t size, const Cancellation &cancellation):
        tree_(size + 1)
    {
        for (std::size_t i = 1; i < tree_.size(); ++i)
        {
            cancellation.tick();
            tree_[i] = static_cast<int>(i & (~i + 1));
        }
    }

    void erase(int position)
    {
        for (std::size_t i = static_cast<std::size_t>(position) + 1;
             i < tree_.size(); i += i & (~i + 1)) --tree_[i];
    }

    int between(int begin, int end) const { return prefix(end) - prefix(begin); }
};

struct Token
{
    Fragment fragment;
    int previous = -1;
    int next = -1;
    int group = -1;
};

struct Run
{
    int last = -1;
    int count = 0;
};

struct Group
{
    int symbol = -1;
    int nonoverlapping = 0;
    std::size_t version = 0;
    // A run consists of consecutive edges of the active-token path.
    std::map<int, Run> runs;
};

struct Candidate
{
    int saving;
    int first;
    int group;
    std::size_t version;

    bool operator<(const Candidate &other) const
    {
        if (saving != other.saving) return saving < other.saving;
        if (first != other.first) return first > other.first;
        return group > other.group;
    }
};

class Compressor
{
    Result result_;
    Cancellation cancellation_;
    Substrings substrings_;
    ActivePositions active_;
    std::vector<Token> tokens_;
    std::map<Key, int> groupIds_;
    std::vector<Group> groups_;
    std::priority_queue<Candidate> candidates_;

    void changed(int id)
    {
        auto &group = groups_[static_cast<std::size_t>(id)];
        ++group.version;
        const int saving = group.nonoverlapping - (group.symbol < 0 ? 1 : 0);
        if (saving > 0)
            candidates_.push({saving, group.runs.begin()->first, id, group.version});
    }

    void removePair(int position)
    {
        if (position < 0) return;
        auto &token = tokens_[static_cast<std::size_t>(position)];
        if (token.group < 0) return;
        const int id = token.group;
        auto &group = groups_[static_cast<std::size_t>(id)];
        auto run = std::prev(group.runs.upper_bound(position));
        const int first = run->first;
        const Run old = run->second;
        const int before = active_.between(first, position);
        const int after = old.count - before - 1;
        group.nonoverlapping -= (old.count + 1) / 2;
        group.runs.erase(run);
        if (before != 0)
        {
            group.runs.emplace(first, Run{token.previous, before});
            group.nonoverlapping += (before + 1) / 2;
        }
        if (after != 0)
        {
            group.runs.emplace(token.next, Run{old.last, after});
            group.nonoverlapping += (after + 1) / 2;
        }
        token.group = -1;
        changed(id);
    }

    void addPair(int position)
    {
        if (position < 0) return;
        auto &token = tokens_[static_cast<std::size_t>(position)];
        if (token.next < 0) return;
        const auto &right = tokens_[static_cast<std::size_t>(token.next)];
        const auto key = substrings_.canonical(token.fragment.offset,
            token.fragment.length + right.fragment.length).first;
        auto entry = groupIds_.lower_bound(key);
        if (entry == groupIds_.end() || key < entry->first)
        {
            if (groups_.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::length_error("too many distinct string Re-Pair products to index");
            entry = groupIds_.emplace_hint(entry, key, static_cast<int>(groups_.size()));
            groups_.emplace_back();
        }
        const int id = entry->second;
        auto &group = groups_[static_cast<std::size_t>(id)];
        int first = position;
        Run merged{position, 1};
        if (token.previous >= 0 &&
            tokens_[static_cast<std::size_t>(token.previous)].group == id)
        {
            auto leftRun = std::prev(group.runs.upper_bound(token.previous));
            first = leftRun->first;
            merged.count += leftRun->second.count;
            group.nonoverlapping -= (leftRun->second.count + 1) / 2;
            group.runs.erase(leftRun);
        }
        if (right.group == id)
        {
            auto rightRun = group.runs.find(token.next);
            merged.last = rightRun->second.last;
            merged.count += rightRun->second.count;
            group.nonoverlapping -= (rightRun->second.count + 1) / 2;
            group.runs.erase(rightRun);
        }
        group.runs.emplace(first, merged);
        group.nonoverlapping += (merged.count + 1) / 2;
        token.group = id;
        changed(id);
    }

    std::vector<int> selectedPairs(int id) const
    {
        const auto &group = groups_[static_cast<std::size_t>(id)];
        std::vector<int> result;
        result.reserve(static_cast<std::size_t>(group.nonoverlapping));
        for (const auto &[first, run] : group.runs)
        {
            int current = first;
            for (int remaining = run.count; remaining > 0; remaining -= 2)
            {
                cancellation_.tick();
                result.push_back(current);
                if (remaining > 2)
                    current = tokens_[static_cast<std::size_t>(
                        tokens_[static_cast<std::size_t>(current)].next)].next;
            }
        }
        return result;
    }

    void defineRule(int id, int representative)
    {
        auto &group = groups_[static_cast<std::size_t>(id)];
        if (group.symbol >= 0) return;
        const auto &first = tokens_[static_cast<std::size_t>(representative)];
        const auto &second = tokens_[static_cast<std::size_t>(first.next)];
        const int length = first.fragment.length + second.fragment.length;
        const bool reverse = substrings_.canonical(first.fragment.offset, length).second;
        const auto &left = reverse ? second.fragment : first.fragment;
        const auto &right = reverse ? first.fragment : second.fragment;
        const auto oriented = [&](const Fragment &fragment) {
            return (fragment.reversed != reverse) && !substrings_.palindrome(fragment);
        };
        const auto symbolCount = result_.terminals.size() + result_.rules.size();
        if (symbolCount >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("too many string Re-Pair symbols to index");
        const int symbol = static_cast<int>(symbolCount);
        group.symbol = symbol;
        result_.rules.push_back({symbol, left.symbol, right.symbol,
            oriented(left), oriented(right), length});
        changed(id);
    }

    void replacePair(int position, int symbol)
    {
        auto &left = tokens_[static_cast<std::size_t>(position)];
        const int rightPosition = left.next;
        auto &right = tokens_[static_cast<std::size_t>(rightPosition)];
        const int before = left.previous;
        const int after = right.next;
        // All old edges touching either token must disappear before the
        // topology and active-position counts change.
        removePair(before);
        removePair(position);
        removePair(rightPosition);
        left.fragment.symbol = symbol;
        left.fragment.length += right.fragment.length;
        left.fragment.reversed = substrings_.canonical(
            left.fragment.offset, left.fragment.length).second;
        left.next = after;
        if (after >= 0) tokens_[static_cast<std::size_t>(after)].previous = position;
        right.previous = right.next = -1;
        active_.erase(rightPosition);
        addPair(before);
        addPair(position);
    }

public:
    Compressor(const std::u32string &input, bool acceptReversed,
               std::function<bool()> cancellationRequested = {}):
        cancellation_(std::move(cancellationRequested)),
        substrings_(input, acceptReversed, cancellation_),
        active_(input.size(), cancellation_)
    {
        result_.acceptReversed = acceptReversed;
        result_.trivialUpperBound = static_cast<int>(input.size()) - 1;
        std::map<char32_t, int> terminals;
        tokens_.reserve(input.size());
        for (std::size_t i = 0; i < input.size(); ++i)
        {
            cancellation_.tick();
            const auto [entry, inserted] = terminals.emplace(input[i],
                static_cast<int>(terminals.size()));
            if (inserted) result_.terminals.push_back({entry->second, input[i]});
            const int position = static_cast<int>(i);
            tokens_.push_back({{entry->second, false, position, 1}, position - 1,
                i + 1 < input.size() ? position + 1 : -1, -1});
        }
        for (std::size_t i = 0; i + 1 < input.size(); ++i)
        {
            cancellation_.tick();
            addPair(static_cast<int>(i));
        }
    }

    Result run()
    {
        cancellation_.poll();
        while (!candidates_.empty())
        {
            cancellation_.tick();
            const auto candidate = candidates_.top();
            candidates_.pop();
            if (groups_[static_cast<std::size_t>(candidate.group)].version !=
                candidate.version) continue;
            const auto pairs = selectedPairs(candidate.group);
            defineRule(candidate.group, pairs.front());
            const int symbol = groups_[static_cast<std::size_t>(candidate.group)].symbol;
            for (const int position : pairs)
            {
                cancellation_.tick();
                replacePair(position, symbol);
            }
        }
        for (int position = tokens_.empty() ? -1 : 0; position >= 0;
             position = tokens_[static_cast<std::size_t>(position)].next)
        {
            cancellation_.tick();
            result_.residual.push_back(tokens_[static_cast<std::size_t>(position)].fragment);
        }
        cancellation_.poll();
        result_.ruleCount = static_cast<int>(result_.rules.size());
        result_.remainingFragments = static_cast<int>(result_.residual.size());
        result_.upperBound = result_.ruleCount + result_.remainingFragments - 1;
        return std::move(result_);
    }
};
} // namespace implementation

/**
 * Calculate a deterministic constructive bound, without Unicode normalization.
 * Malformed UTF-8 and strings too large for scalar indices throw invalid_argument.
 */
inline Result calculate(std::string_view input, bool acceptReversed = false)
{
    const auto decoded = stringEncoding::decodeInput(input);
    return implementation::Compressor(decoded, acceptReversed).run();
}

/**
 * Write an ASCII JSON certificate for a result belonging to input. Expanding
 * terminals and topologically ordered rules, orienting residual references,
 * and concatenating them reconstructs input without consulting this algorithm.
 */
inline void writeJson(const Result &result, std::string_view input, std::ostream &output)
{
    const auto decoded = stringEncoding::decodeInput(input);
    output << "{\"schema\":\"string-repair-assembly-v1\",\"upper_bound\":" << result.upperBound
           << ",\"trivial_upper_bound\":" << result.trivialUpperBound
           << ",\"rule_count\":" << result.ruleCount
           << ",\"remaining_fragments\":" << result.remainingFragments
           << ",\"accept_reversed\":" << (result.acceptReversed ? "true" : "false")
           << ",\"length\":" << decoded.size() << ",\"input\":";
    stringEncoding::writeJsonString(input, output);
    output << ",\"terminals\":[";
    for (std::size_t i = 0; i < result.terminals.size(); ++i)
    {
        if (i != 0) output << ',';
        const auto &terminal = result.terminals[i];
        output << "{\"id\":" << terminal.id << ",\"code_point\":"
               << static_cast<std::uint32_t>(terminal.codePoint) << '}';
    }
    output << "],\"rules\":[";
    for (std::size_t i = 0; i < result.rules.size(); ++i)
    {
        if (i != 0) output << ',';
        const auto &rule = result.rules[i];
        output << "{\"id\":" << rule.id << ",\"left\":" << rule.left
               << ",\"right\":" << rule.right
               << ",\"left_reversed\":" << (rule.leftReversed ? "true" : "false")
               << ",\"right_reversed\":" << (rule.rightReversed ? "true" : "false")
               << ",\"length\":" << rule.length << '}';
    }
    output << "],\"residual\":[";
    for (std::size_t i = 0; i < result.residual.size(); ++i)
    {
        if (i != 0) output << ',';
        const auto &fragment = result.residual[i];
        output << "{\"symbol\":" << fragment.symbol
               << ",\"reversed\":" << (fragment.reversed ? "true" : "false")
               << ",\"offset\":" << fragment.offset << ",\"length\":" << fragment.length << '}';
    }
    output << "]}";
}
} // namespace parallelassemblycpp::detail::stringRepair
