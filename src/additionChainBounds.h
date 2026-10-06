#ifndef PARALLELASSEMBLYCPP_ADDITION_CHAIN_BOUNDS_H
#define PARALLELASSEMBLYCPP_ADDITION_CHAIN_BOUNDS_H

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace assembly_bounds
{

#ifdef PARALLELASSEMBLYCPP_DISABLE_VECTOR_CHAIN_BOUNDS
inline constexpr bool vectorBoundsEnabled = false;
#else
inline constexpr bool vectorBoundsEnabled = true;
#endif

namespace detail
{

inline int ceilLog2(std::uint64_t value) noexcept
{
    if (value <= 1) return 0;
    return std::bit_width(value - 1);
}

// Exact unrestricted scalar addition-chain lengths, indexed by target.
// The zero entry is the empty-object convention. These constants are verified
// exhaustively (all earlier-pair sums, not just star chains) by the unit test.
// Keeping the table immutable avoids startup search and worker-local copies.
inline constexpr std::array<unsigned char, 257> scalarLengths = {
    0, 0, 1, 2, 2, 3, 3, 4, 3, 4, 4, 5, 4, 5, 5, 5,
    4, 5, 5, 6, 5, 6, 6, 6, 5, 6, 6, 6, 6, 7, 6, 7,
    5, 6, 6, 7, 6, 7, 7, 7, 6, 7, 7, 7, 7, 7, 7, 8,
    6, 7, 7, 7, 7, 8, 7, 8, 7, 8, 8, 8, 7, 8, 8, 8,
    6, 7, 7, 8, 7, 8, 8, 9, 7, 8, 8, 8, 8, 8, 8, 9,
    7, 8, 8, 8, 8, 8, 8, 9, 8, 9, 8, 9, 8, 9, 9, 9,
    7, 8, 8, 8, 8, 9, 8, 9, 8, 9, 9, 9, 8, 9, 9, 9,
    8, 9, 9, 9, 9, 9, 9, 9, 8, 9, 9, 9, 9, 9, 9, 10,
    7, 8, 8, 9, 8, 9, 9, 9, 8, 9, 9, 10, 9, 10, 10, 10,
    8, 9, 9, 9, 9, 9, 9, 10, 9, 9, 9, 10, 9, 10, 10, 10,
    8, 9, 9, 9, 9, 9, 9, 10, 9, 10, 9, 10, 9, 10, 10, 10,
    9, 10, 10, 10, 9, 10, 10, 10, 9, 10, 10, 10, 10, 10, 10, 11,
    8, 9, 9, 9, 9, 10, 9, 10, 9, 10, 10, 10, 9, 10, 10, 10,
    9, 10, 10, 10, 10, 10, 10, 10, 9, 10, 10, 10, 10, 10, 10, 11,
    9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 10, 11, 10, 11,
    9, 10, 10, 10, 10, 10, 10, 11, 10, 10, 10, 11, 10, 11, 11, 10,
    8
};

inline int scalarLowerBound64(std::uint64_t value) noexcept
{
    return value < scalarLengths.size()
        ? scalarLengths[static_cast<std::size_t>(value)] : ceilLog2(value);
}

} // namespace detail

/** A certified lower bound; exact for targets 0 through 256. */
inline int scalarLowerBound(int target) noexcept
{
    return target <= 1 ? 0 : detail::scalarLowerBound64(
        static_cast<std::uint64_t>(target)
    );
}

/**
 * Bounds for vector addition chains starting from the nonzero unit vectors.
 *
 * Coordinate names are irrelevant, so permutations and inserted zeros share a
 * cache entry. Summing the coordinates maps any vector chain to a scalar chain
 * after redundant values are removed. Connecting k different primitive types
 * requires at least k-1 additions, or k when any coordinate repeats. These
 * are lower bounds, not lengths of heuristic construction paths.
 *
 * Small targets additionally use exhaustive bounded search. A depth contributes
 * a stronger bound only after it has been completely disproved; exhausting the
 * work budget never implies that the depth was impossible. The budget counts
 * every pair considered, including rejected pairs, across all searched depths.
 *
 * Instances own all mutable state: use one instance per search/worker, or lock
 * externally when sharing. Cache eviction is deterministic and capacity bounded.
 */
class VectorBoundCache
{
public:
    static constexpr int maximumExactTotal = 16;
    static constexpr std::size_t defaultCapacity = 256;
    static constexpr std::size_t defaultWorkLimit = 4000;

    explicit VectorBoundCache(
        std::size_t capacity = defaultCapacity,
        std::size_t workLimit = defaultWorkLimit
    ): capacity_(capacity), workLimit_(workLimit)
    {}

    /**
     * Stop exact-search refinement once a certificate reaches sufficientBound.
     * Cheap projections may already exceed that threshold. Cached certificates
     * remain valid across thresholds, but an earlier stopped query can leave a
     * weaker certificate than a fresh full-budget query would produce.
     */
    int lowerBound(
        std::vector<int> counts,
        int sufficientBound = std::numeric_limits<int>::max()
    )
    {
        std::uint64_t total = 0;
        for (const int count : counts)
        {
            if (count < 0)
                throw std::invalid_argument("negative addition-chain count");
            total += static_cast<std::uint64_t>(count);
        }
        counts.erase(std::remove(counts.begin(), counts.end(), 0), counts.end());
        if (counts.empty()) return 0;
        if (counts.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("too many addition-chain coordinates");
        std::sort(counts.begin(), counts.end());
        const auto cached = cache_.find(counts);
        if (cached != cache_.end()) return cached->second;

        // With k unit-vector sources and k-1 binary additions, the ancestor
        // graph has exactly enough edges to be a tree. There can be no shared
        // operand, so each source contributes once. A repeated coordinate thus
        // needs at least one further addition. Counts are sorted and nonzero.
        const int supportBound = static_cast<int>(counts.size()) -
            (counts.back() == 1 ? 1 : 0);
        int bound = std::max(detail::scalarLowerBound64(total), supportBound);
        // Projection onto a single coordinate is also a scalar chain once
        // zeros and repeated values are removed. Exact scalar lengths are not
        // monotone, so this can improve on the total-count projection.
        for (const int count : counts)
            bound = std::max(bound, scalarLowerBound(count));

        if (counts.size() > 1 && total <= maximumExactTotal
            && bound < static_cast<int>(total) - 1 && workLimit_ != 0
            && bound < sufficientBound)
        {
            ExactSearch search(counts, workLimit_);
            while (bound < static_cast<int>(total) - 1 && bound < sufficientBound)
            {
                const SearchResult result = search.run(bound);
                if (result != SearchResult::impossible) break;
                ++bound;
            }
        }

        if (capacity_ != 0)
        {
            if (cache_.size() == capacity_) cache_.erase(cache_.begin());
            cache_.emplace(std::move(counts), bound);
        }
        return bound;
    }

    [[nodiscard]] std::size_t size() const noexcept { return cache_.size(); }
    void clear() noexcept { cache_.clear(); }

private:
    enum class SearchResult { found, impossible, aborted };

    struct Element
    {
        std::array<unsigned char, maximumExactTotal> counts{};
        std::uint32_t code = 0;
        unsigned int weight = 0;
    };

    class ExactSearch
    {
    public:
        ExactSearch(const std::vector<int> &counts, std::size_t workLimit):
            dimensions_(counts.size()), workLeft_(workLimit)
        {
            std::uint32_t multiplier = 1;
            for (std::size_t i = 0; i < dimensions_ && i < maximumExactTotal; ++i)
            {
                target_.counts[i] = static_cast<unsigned char>(counts[i]);
                target_.weight += static_cast<unsigned int>(counts[i]);
                target_.code += static_cast<std::uint32_t>(counts[i]) * multiplier;
                Element primitive;
                primitive.counts[i] = 1;
                primitive.code = multiplier;
                primitive.weight = 1;
                chain_.push_back(primitive);
                multiplier *= static_cast<std::uint32_t>(counts[i] + 1);
            }
            chain_.reserve(dimensions_ + maximumExactTotal);
        }

        SearchResult run(int additions)
        {
            return visit(additions, 0);
        }

    private:
        SearchResult visit(int remaining, std::uint32_t lastGenerated)
        {
            if (remaining == 0) return SearchResult::impossible;
            unsigned int maximumWeight = 1;
            std::array<unsigned char, maximumExactTotal> maximumCounts{};
            for (const Element &value : chain_)
            {
                maximumWeight = std::max(maximumWeight, value.weight);
                for (std::size_t c = 0; c < dimensions_ && c < maximumExactTotal; ++c)
                    maximumCounts[c] = std::max(maximumCounts[c], value.counts[c]);
            }
            // Each addition can at most double any existing coordinate or the
            // total weight. These tests only discard provably unreachable goals.
            if ((maximumWeight << remaining) < target_.weight)
                return SearchResult::impossible;
            for (std::size_t c = 0; c < dimensions_ && c < maximumExactTotal; ++c)
                if ((static_cast<unsigned int>(maximumCounts[c]) << remaining)
                    < target_.counts[c])
                    return SearchResult::impossible;

            std::vector<Element> candidates;
            candidates.reserve(chain_.size() * (chain_.size() + 1) / 2);
            for (std::size_t i = 0; i < chain_.size(); ++i)
                for (std::size_t j = 0; j <= i; ++j)
                {
                    if (workLeft_ == 0) return SearchResult::aborted;
                    --workLeft_;
                    Element candidate;
                    candidate.code = chain_[i].code + chain_[j].code;
                    candidate.weight = chain_[i].weight + chain_[j].weight;
                    if (candidate.code <= lastGenerated
                        || candidate.code > target_.code
                        || (remaining == 1 && candidate.code != target_.code))
                        continue;
                    bool fits = true;
                    for (std::size_t c = 0; c < dimensions_ && c < maximumExactTotal; ++c)
                    {
                        candidate.counts[c] = static_cast<unsigned char>(
                            chain_[i].counts[c] + chain_[j].counts[c]
                        );
                        if (candidate.counts[c] > target_.counts[c])
                        {
                            fits = false;
                            break;
                        }
                    }
                    if (!fits) continue;
                    if (candidate.code == target_.code) return SearchResult::found;
                    if (remaining > 1) candidates.push_back(candidate);
                }
            std::sort(candidates.begin(), candidates.end(),
                [](const Element &left, const Element &right)
                { return left.code > right.code; });
            candidates.erase(std::unique(candidates.begin(), candidates.end(),
                [](const Element &left, const Element &right)
                { return left.code == right.code; }), candidates.end());
            for (const Element &candidate : candidates)
            {
                chain_.push_back(candidate);
                const SearchResult result = visit(remaining - 1, candidate.code);
                chain_.pop_back();
                if (result != SearchResult::impossible) return result;
            }
            return SearchResult::impossible;
        }

        // Mixed-radix encoding is injective for componentwise valid vectors and
        // additive for sums that fit the target. Its maximum is at most 2^16-1.
        // Every nonprimitive element can be ordered by increasing code because
        // both operands have strictly smaller positive codes. Enumerating that
        // order therefore covers every chain, including non-star chains.
        std::size_t dimensions_;
        std::size_t workLeft_;
        Element target_;
        std::vector<Element> chain_;
    };

    std::size_t capacity_;
    std::size_t workLimit_;
    std::map<std::vector<int>, int> cache_;
};

} // namespace assembly_bounds

#endif
