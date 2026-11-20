#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace parallelassemblycpp::detail
{
// Non-overlapping LZ factorisation as used by assemblycfg (MIT, see
// docs/assemblycfg-license.md). A left-to-right traversal of a binary grammar
// expands each rule only on first use. Its other occurrences are copies from
// the completed prefix. Thus a grammar with k joins gives at most k+1 phrases.
// Greedy longest phrases minimize this count: starting farther along can
// always copy the remaining suffix of any phrase reachable from an earlier cut.
template<typename Tick>
int directedLzLowerBound(std::u32string_view input, Tick tick)
{
    if (input.size() < 2) return static_cast<int>(input.size()) - 1;
    constexpr auto missing = std::numeric_limits<std::size_t>::max();
    struct State
    {
        std::size_t length = 0;
        std::size_t link = std::numeric_limits<std::size_t>::max();
        std::unordered_map<char32_t, std::size_t> edges;
    };
    std::vector<State> states(1);
    std::size_t last = 0;
    // Extend the suffix automaton of the already consumed prefix by one scalar.
    const auto append = [&](char32_t symbol) {
        const auto current = states.size();
        states.push_back({states[last].length + 1, 0, {}});
        auto parent = last;
        while (parent != missing && states[parent].edges.emplace(symbol, current).second)
        {
            tick();
            parent = states[parent].link;
        }
        if (parent != missing)
        {
            const auto child = states[parent].edges.at(symbol);
            if (states[parent].length + 1 == states[child].length)
                states[current].link = child;
            else
            {
                const auto clone = states.size();
                states.push_back(states[child]);
                states[clone].length = states[parent].length + 1;
                states[child].link = states[current].link = clone;
                while (parent != missing)
                {
                    tick();
                    auto &edge = states[parent].edges.at(symbol);
                    if (edge != child) break;
                    edge = clone;
                    parent = states[parent].link;
                }
            }
        }
        last = current;
    };

    int joins = -1;
    for (std::size_t start = 0; start < input.size(); ++joins)
    {
        std::size_t state = 0;
        auto end = start;
        while (end < input.size())
        {
            tick();
            const auto next = states[state].edges.find(input[end]);
            if (next == states[state].edges.end()) break;
            state = next->second;
            ++end;
        }
        end = std::max(end, start + 1);
        while (start < end)
        {
            tick();
            append(input[start++]);
        }
    }
    return joins;
}
}
