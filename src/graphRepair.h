#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <ostream>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * GraphRePair-inspired molecular assembly upper bound.
 *
 * Active tokens partition the physical bonds into connected fragments. A
 * production joins two incident, edge-disjoint tokens. Exact canonical forms
 * of the FULL expanded labelled graph identify reusable products, including
 * joins at multiple vertices. All vertices remain available as attachment
 * sites; this deliberately relaxes GraphRePair's fixed boundary interfaces to
 * match molecular assembly. It is not the paper's linear-time compressor.
 *
 * A new binary production costs one join; using an existing production costs
 * zero. Replacing k disjoint pairs therefore saves k-1 (or k for an existing
 * product) against the residual join count. At termination, r connected tokens
 * covering c nonempty components can be joined in r-c steps. The repository's
 * default disconnected convention adds c-1, giving rules+r-1. Optional
 * compensation omits those c-1 steps. Isolated atoms are not bond primitives.
 *
 * Canonicalisation uses the solver's thread-local interner. Like the public
 * calculation API, callers must not overlap calculations in one thread.
 */
namespace graphRepair
{
struct Fragment
{
    int symbol = -1;
    std::vector<int> edges;
};

struct Rule
{
    int id = -1;
    int left = -1;
    int right = -1;
    std::vector<int> leftEdges;
    std::vector<int> rightEdges;
    std::vector<int> edges;
};

struct Result
{
    int upperBound = 0;
    int trivialUpperBound = 0;
    int ruleCount = 0;
    int remainingFragments = 0;
    int components = 0;
    bool compensateDisjoint = false;
    std::vector<Fragment> terminals;
    std::vector<Rule> rules;
    std::vector<Fragment> residual;
};

namespace implementation
{
#ifdef PARALLELASSEMBLYCPP_LIBRARY_BUILD
using CanonicalMap = std::unordered_map<graphHash, int, graphHashHasher>;
#else
using CanonicalMap = std::unordered_map<graphHash, int>;
#endif

struct Token
{
    Fragment fragment;
    std::vector<int> vertices;
};

struct Pair
{
    std::size_t first;
    std::size_t second;
};

struct Group
{
    int symbol = -1;
    std::vector<Pair> occurrences;
};

inline std::vector<int> merge(const std::vector<int> &left,
                              const std::vector<int> &right)
{
    std::vector<int> result;
    result.reserve(left.size() + right.size());
    std::set_union(left.begin(), left.end(), right.begin(), right.end(),
                   std::back_inserter(result));
    return result;
}

inline bool incident(const std::vector<int> &left,
                     const std::vector<int> &right)
{
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < left.size() && j < right.size())
    {
        if (left[i] == right[j]) return true;
        if (left[i] < right[j]) ++i;
        else ++j;
    }
    return false;
}

/** Local graph construction retains every atom and bond label exactly. */
inline graphHash canonical(const molGraph &source,
                           const std::vector<MoleculeEdge> &edges,
                           const std::vector<int> &selected)
{
    molGraph fragment;
    std::vector<int> local(source.atoms.size(), -1);
    for (int edgeId : selected)
    {
        const auto &edge = edges[static_cast<std::size_t>(edgeId)];
        for (int vertex : {edge.sourceAtomIndex, edge.targetAtomIndex})
        {
            if (local[static_cast<std::size_t>(vertex)] != -1) continue;
            local[static_cast<std::size_t>(vertex)] =
                static_cast<int>(fragment.atoms.size());
            fragment.addAtom(source.atoms[static_cast<std::size_t>(vertex)].atomType);
        }
        fragment.addBond(local[static_cast<std::size_t>(edge.sourceAtomIndex)],
                         local[static_cast<std::size_t>(edge.targetAtomIndex)],
                         source.bondType(edge.sourceAtomIndex, edge.sourceBondIndex));
    }
    const bool cyclic = selected.size() >= fragment.atoms.size();
    return graphHash(fragment, cyclic);
}

inline int componentCount(const molGraph &graph)
{
    std::vector<bool> visited(graph.atoms.size(), false);
    std::vector<std::size_t> pending;
    int components = 0;
    for (std::size_t vertex = 0; vertex < graph.atoms.size(); ++vertex)
    {
        if (visited[vertex] || graph.atoms[vertex].bonds.empty()) continue;
        ++components;
        visited[vertex] = true;
        pending.push_back(vertex);
        while (!pending.empty())
        {
            const auto current = pending.back();
            pending.pop_back();
            for (const bond &edge : graph.atoms[current].bonds)
            {
                const auto neighbour = static_cast<std::size_t>(edge.neighbourAtomIndex);
                if (visited[neighbour]) continue;
                visited[neighbour] = true;
                pending.push_back(neighbour);
            }
        }
    }
    return components;
}
} // namespace implementation

/** Calculate a constructive bound without enumerating all connected subgraphs. */
inline Result calculate(const molGraph &graph, bool compensateDisjoint = false)
{
    using namespace implementation;
    Result result;
    result.compensateDisjoint = compensateDisjoint;
    result.components = componentCount(graph);
    const auto edges = graph.writeEdgeList();
    if (edges.empty()) return result;
    const int adjustment = compensateDisjoint ? result.components : 1;
    result.trivialUpperBound = static_cast<int>(edges.size()) - adjustment;

    // These keys never outlive this invocation or its canonical interner.
    CanonicalMap dictionary;
    std::vector<Token> tokens;
    std::vector<std::size_t> active;
    tokens.reserve(edges.size() * 2);
    active.reserve(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        const std::vector<int> selected{static_cast<int>(i)};
        auto key = canonical(graph, edges, selected);
        const auto [position, inserted] = dictionary.emplace(
            std::move(key), static_cast<int>(dictionary.size()));
        if (inserted) result.terminals.push_back({position->second, selected});
        std::vector<int> vertices{edges[i].sourceAtomIndex, edges[i].targetAtomIndex};
        std::sort(vertices.begin(), vertices.end());
        tokens.push_back({{position->second, selected}, std::move(vertices)});
        active.push_back(i);
    }

    // Stable token IDs let unchanged adjacent pairs reuse their exact forms.
    std::map<std::pair<std::size_t, std::size_t>, graphHash> pairKeys;
    while (active.size() > 1)
    {
        CanonicalMap groupIds;
        std::vector<Group> groups;
        for (std::size_t i = 0; i < active.size(); ++i)
        {
            for (std::size_t j = i + 1; j < active.size(); ++j)
            {
                const auto first = active[i];
                const auto second = active[j];
                if (!incident(tokens[first].vertices, tokens[second].vertices)) continue;
                const auto pairId = std::minmax(first, second);
                auto cached = pairKeys.find(pairId);
                if (cached == pairKeys.end())
                {
                    auto selected = merge(tokens[first].fragment.edges,
                                          tokens[second].fragment.edges);
                    cached = pairKeys.emplace(pairId,
                        canonical(graph, edges, selected)).first;
                }
                const auto [group, inserted] = groupIds.emplace(
                    cached->second, static_cast<int>(groups.size()));
                if (inserted)
                {
                    const auto known = dictionary.find(cached->second);
                    groups.push_back({known == dictionary.end() ? -1 : known->second, {}});
                }
                groups[static_cast<std::size_t>(group->second)].occurrences.push_back(
                    {first, second});
            }
        }

        // Greedy edge-disjoint occurrence selection. Overlap at vertices is
        // permitted; overlap at tokens (and hence physical bonds) is not.
        int bestSaving = 0;
        int bestGroup = -1;
        std::vector<Pair> bestPairs;
        for (std::size_t groupId = 0; groupId < groups.size(); ++groupId)
        {
            std::vector<bool> used(tokens.size(), false);
            std::vector<Pair> selected;
            for (const auto &occurrence : groups[groupId].occurrences)
            {
                if (used[occurrence.first] || used[occurrence.second]) continue;
                used[occurrence.first] = used[occurrence.second] = true;
                selected.push_back(occurrence);
            }
            const int saving = static_cast<int>(selected.size()) -
                (groups[groupId].symbol == -1 ? 1 : 0);
            if (saving > bestSaving)
            {
                bestSaving = saving;
                bestGroup = static_cast<int>(groupId);
                bestPairs = std::move(selected);
            }
        }
        if (bestGroup == -1) break;
        auto &group = groups[static_cast<std::size_t>(bestGroup)];
        if (group.symbol == -1)
        {
            const auto representative = bestPairs.front();
            const auto &left = tokens[representative.first].fragment;
            const auto &right = tokens[representative.second].fragment;
            group.symbol = static_cast<int>(dictionary.size());
            dictionary.emplace(pairKeys.at(std::minmax(representative.first,
                                                       representative.second)), group.symbol);
            result.rules.push_back({group.symbol, left.symbol, right.symbol,
                                    left.edges, right.edges, merge(left.edges, right.edges)});
        }

        std::vector<bool> removed(tokens.size(), false);
        std::vector<std::size_t> next;
        next.reserve(active.size() - bestPairs.size());
        for (const auto &occurrence : bestPairs)
        {
            const auto &left = tokens[occurrence.first];
            const auto &right = tokens[occurrence.second];
            Token replacement{{group.symbol, merge(left.fragment.edges, right.fragment.edges)},
                              merge(left.vertices, right.vertices)};
            removed[occurrence.first] = removed[occurrence.second] = true;
            next.push_back(tokens.size());
            tokens.push_back(std::move(replacement));
        }
        for (auto token : active) if (!removed[token]) next.push_back(token);
        std::sort(next.begin(), next.end(), [&](std::size_t left, std::size_t right) {
            return tokens[left].fragment.edges.front() < tokens[right].fragment.edges.front();
        });
        active = std::move(next);
    }
    for (auto token : active) result.residual.push_back(std::move(tokens[token].fragment));
    result.ruleCount = static_cast<int>(result.rules.size());
    result.remainingFragments = static_cast<int>(result.residual.size());
    result.upperBound = result.ruleCount + result.remainingFragments - adjustment;
    return result;
}

/** A replayable certificate: edge IDs refer to the listed original graph. */
inline void writeJson(const Result &result, const molGraph &graph, std::ostream &output)
{
    const auto list = [&](const std::vector<int> &values) {
        output << '[';
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            if (i != 0) output << ',';
            output << values[i];
        }
        output << ']';
    };
    const auto fragments = [&](const std::vector<Fragment> &values) {
        output << '[';
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            if (i != 0) output << ',';
            output << "{\"" << (&values == &result.terminals ? "id" : "symbol")
                   << "\":" << values[i].symbol << ",\"edges\":";
            list(values[i].edges);
            output << '}';
        }
        output << ']';
    };
    output << "{\"schema\":\"graph-repair-assembly-v1\",\"upper_bound\":" << result.upperBound
           << ",\"trivial_upper_bound\":" << result.trivialUpperBound
           << ",\"rule_count\":" << result.ruleCount
           << ",\"remaining_fragments\":" << result.remainingFragments
           << ",\"components\":" << result.components
           << ",\"compensate_disjoint\":" << (result.compensateDisjoint ? "true" : "false")
           << ",\"atoms\":[";
    for (std::size_t i = 0; i < graph.atoms.size(); ++i)
    {
        if (i != 0) output << ',';
        parallelassemblycpp::detail::stringAssembly::implementation::writeJsonString(
            graph.atoms[i].atomType, output);
    }
    output << "],\"edges\":[";
    const auto edges = graph.writeEdgeList();
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        if (i != 0) output << ',';
        const auto &edge = edges[i];
        output << '[' << edge.sourceAtomIndex << ',' << edge.targetAtomIndex << ','
               << graph.bondType(edge.sourceAtomIndex, edge.sourceBondIndex) << ']';
    }
    output << "],\"terminals\":";
    fragments(result.terminals);
    output << ",\"rules\":[";
    for (std::size_t i = 0; i < result.rules.size(); ++i)
    {
        if (i != 0) output << ',';
        const auto &rule = result.rules[i];
        output << "{\"id\":" << rule.id << ",\"left\":" << rule.left
               << ",\"right\":" << rule.right << ",\"left_edges\":";
        list(rule.leftEdges);
        output << ",\"right_edges\":";
        list(rule.rightEdges);
        output << ",\"edges\":";
        list(rule.edges);
        output << '}';
    }
    output << "],\"residual\":";
    fragments(result.residual);
    output << '}';
}
} // namespace graphRepair
