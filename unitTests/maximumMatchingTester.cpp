#include "../src/maximumMatching.h"

#include <bit>
#include <iostream>
#include <random>
#include <stdexcept>

namespace
{
using Graph = std::vector<std::vector<std::size_t>>;

int optimum(const Graph &neighbors, unsigned int available, std::vector<int> &memo)
{
    if (available == 0) return 0;
    int &best = memo[available];
    if (best >= 0) return best;
    const auto first = static_cast<std::size_t>(std::countr_zero(available));
    const auto rest = available & (available - 1);
    best = optimum(neighbors, rest, memo);
    for (const auto next : neighbors[first])
    {
        const auto bit = 1U << next;
        if (rest & bit) best = std::max(best, 1 + optimum(neighbors, rest ^ bit, memo));
    }
    return best;
}

void check(const Graph &neighbors)
{
    const auto size = neighbors.size();
    std::vector<std::size_t> greedy(size, size);
    for (std::size_t i = 0; i < size; ++i)
        for (auto j : neighbors[i])
            if (greedy[i] == size && greedy[j] == size)
            {
                greedy[i] = j;
                greedy[j] = i;
            }
    const auto matching = assembly_matching::MaximumMatching(neighbors, std::move(greedy)).run();
    int matched = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
        const auto j = matching[i];
        if (j == size) continue;
        if (matching.at(j) != i ||
            std::find(neighbors[i].begin(), neighbors[i].end(), j) == neighbors[i].end())
            throw std::logic_error("invalid matching");
        ++matched;
    }
    std::vector<int> memo(1U << size, -1);
    if (matched != 2 * optimum(neighbors, (1U << size) - 1, memo))
        throw std::logic_error("matching is too small");
}
}

int main()
{
    // All simple graphs on six vertices cover alternating paths, odd cycles
    // and disconnected products. Larger cases exercise nested blossoms.
    for (unsigned int mask = 0; mask < (1U << 15); ++mask)
    {
        Graph neighbors(6);
        unsigned int bit = 1;
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = i + 1; j < 6; ++j, bit <<= 1)
                if (mask & bit)
                {
                    neighbors[i].push_back(j);
                    neighbors[j].push_back(i);
                }
        check(neighbors);
    }
    std::mt19937 random(17);
    for (int trial = 0; trial < 150; ++trial)
    {
        Graph neighbors(13);
        for (std::size_t i = 0; i < neighbors.size(); ++i)
            for (std::size_t j = 0; j < i; ++j)
                if (random() % 4 == 0)
                {
                    neighbors[i].push_back(j);
                    neighbors[j].push_back(i);
                }
        check(neighbors);
    }
    std::cout << "Maximum matching tests passed\n";
}
