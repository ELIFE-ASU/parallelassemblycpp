#pragma once

#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#include "../src/molGraph.h"

namespace graphTestFixtures
{

using edgeSpec = std::tuple<int, int, short>;

inline void require(bool condition, const char *message)
{
    if (condition) return;
    std::cerr << "invalid graph fixture: " << message << '\n';
    std::abort();
}

/** Build a labelled graph, preserving the requested edge insertion order. */
inline molGraph makeGraph(
    const std::vector<std::string> &labels,
    const std::vector<edgeSpec> &edges,
    const std::vector<int> &oldToNew = {},
    bool reverseEdges = false
)
{
    require(
        labels.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "graph is too large for integer vertex IDs"
    );
    std::vector<int> permutation = oldToNew;
    if (permutation.empty())
    {
        permutation.resize(labels.size());
        std::iota(permutation.begin(), permutation.end(), 0);
    }
    require(permutation.size() == labels.size(), "invalid permutation size");

    std::vector<unsigned char> seen(labels.size(), 0);
    std::vector<std::string> permutedLabels(labels.size());
    for (std::size_t oldIndex = 0; oldIndex < labels.size(); oldIndex++)
    {
        const int replacement = permutation[oldIndex];
        const std::size_t replacementIndex = static_cast<std::size_t>(replacement);
        require(
            replacement >= 0 &&
                replacementIndex < labels.size() &&
                !seen[replacementIndex],
            "permutation is not a bijection"
        );
        seen[replacementIndex] = 1;
        permutedLabels[replacementIndex] = labels[oldIndex];
    }

    molGraph result;
    for (std::string &label : permutedLabels) result.addAtom(label);
    const auto addOne = [&](const edgeSpec &edge)
    {
        const auto [left, right, bondType] = edge;
        const std::size_t leftIndex = static_cast<std::size_t>(left);
        const std::size_t rightIndex = static_cast<std::size_t>(right);
        require(
            left >= 0 && right >= 0 &&
                leftIndex < permutation.size() &&
                rightIndex < permutation.size(),
            "edge endpoint is outside the graph"
        );
        result.addBond(permutation[leftIndex], permutation[rightIndex], bondType);
    };
    if (reverseEdges)
    {
        for (auto edge = edges.rbegin(); edge != edges.rend(); ++edge)
            addOne(*edge);
    }
    else
    {
        for (const edgeSpec &edge : edges) addOne(edge);
    }
    return result;
}

} // namespace graphTestFixtures
