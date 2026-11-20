#pragma once

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

namespace assembly_matching
{
// Edmonds' alternating-tree algorithm for maximum cardinality matching in
// an undirected simple graph. A same-tree edge between outer vertices closes
// an odd cycle. Treat the cycle as one outer vertex while searching, keeping
// the parent links needed to expand an augmenting path through it.
class MaximumMatching
{
    using Index = std::size_t;
    const std::vector<std::vector<Index>> &neighbors_;
    const Index none_;
    std::vector<Index> mate_, parent_, base_, queue_;
    std::vector<bool> outer_, ancestors_, blossom_;

    Index commonBase(Index left, Index right)
    {
        std::fill(ancestors_.begin(), ancestors_.end(), false);
        for (;;)
        {
            left = base_[left];
            ancestors_[left] = true;
            if (mate_[left] == none_) break;
            left = parent_[mate_[left]];
        }
        right = base_[right];
        while (!ancestors_[right]) right = base_[parent_[mate_[right]]];
        return right;
    }

    void markBranch(Index vertex, Index base, Index child)
    {
        while (base_[vertex] != base)
        {
            blossom_[base_[vertex]] = blossom_[base_[mate_[vertex]]] = true;
            parent_[vertex] = child;
            child = mate_[vertex];
            vertex = parent_[child];
        }
    }

    void enqueue(Index vertex)
    {
        if (outer_[vertex]) return;
        outer_[vertex] = true;
        queue_.push_back(vertex);
    }

    void augmentFrom(Index root)
    {
        std::iota(base_.begin(), base_.end(), Index{0});
        std::fill(parent_.begin(), parent_.end(), none_);
        std::fill(outer_.begin(), outer_.end(), false);
        queue_.clear();
        enqueue(root);
        for (Index head = 0; head < queue_.size(); ++head)
        {
            const auto vertex = queue_[head];
            for (const auto next : neighbors_[vertex])
            {
                if (base_[vertex] == base_[next] || mate_[vertex] == next) continue;
                if (outer_[next])
                {
                    const auto base = commonBase(vertex, next);
                    std::fill(blossom_.begin(), blossom_.end(), false);
                    markBranch(vertex, base, next);
                    markBranch(next, base, vertex);
                    for (Index i = 0; i < none_; ++i)
                        if (blossom_[base_[i]])
                        {
                            base_[i] = base;
                            enqueue(i);
                        }
                }
                else if (parent_[next] == none_)
                {
                    parent_[next] = vertex;
                    if (mate_[next] != none_) enqueue(mate_[next]);
                    else
                    {
                        for (auto end = next; end != none_; )
                        {
                            const auto previous = parent_[end];
                            const auto rest = mate_[previous];
                            mate_[end] = previous;
                            mate_[previous] = end;
                            end = rest;
                        }
                        return;
                    }
                }
            }
        }
    }

public:
    MaximumMatching(const std::vector<std::vector<Index>> &neighbors,
                    std::vector<Index> greedy):
        neighbors_(neighbors), none_(neighbors.size()), mate_(std::move(greedy)),
        parent_(none_), base_(none_), outer_(none_), ancestors_(none_), blossom_(none_)
    {}

    std::vector<Index> run()
    {
        for (Index root = 0; root < none_; ++root)
            if (mate_[root] == none_) augmentFrom(root);
        return std::move(mate_);
    }
};
}
