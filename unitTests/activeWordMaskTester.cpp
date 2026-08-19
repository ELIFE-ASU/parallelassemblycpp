#include "../v5/activeWordMask.h"

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template<typename Left, typename Right>
concept EqualityComparableWith = requires(const Left &left, const Right &right) {
    left == right;
};

template<typename Left, typename Right>
concept BitwiseAndableWith = requires(const Left &left, const Right &right) {
    left & right;
};

static_assert(!std::same_as<EdgeMask, AtomMask>);
static_assert(!std::convertible_to<EdgeMask, AtomMask>);
static_assert(!std::convertible_to<AtomMask, EdgeMask>);
static_assert(!std::assignable_from<EdgeMask &, AtomMask>);
static_assert(!std::assignable_from<AtomMask &, EdgeMask>);
static_assert(!EqualityComparableWith<EdgeMask, AtomMask>);
static_assert(!BitwiseAndableWith<EdgeMask, AtomMask>);

template<typename Mask>
void testWidth(std::size_t width)
{
    Mask::configure(width);
    assert(Mask::size() == width);
    assert(Mask::activeWordCount() == (width + 63) / 64);

    Mask empty;
    assert(empty.none());
    assert(empty.count() == 0);
    assert(empty == 0);
    assert(empty.findFirst() == width);

    Mask full;
    full.set();
    assert(full.count() == width);
    assert(full.all());
    assert((~full).none());

    constexpr std::array<std::size_t, 17> boundaryBits{
        0, 1, 63, 64, 127, 128, 191, 192, 255, 256, 319, 320, 383, 384,
        447, 448, 511
    };
    Mask selected;
    std::vector<std::size_t> expected;
    for (const std::size_t bit : boundaryBits)
    {
        if (bit >= width) continue;
        selected.set(bit);
        expected.push_back(bit);
    }
    assert(selected.count() == expected.size());

    for (
        std::size_t wordIndex = 0;
        wordIndex < Mask::activeWordCount();
        wordIndex++
    )
    {
        unsigned long long expectedWord = 0;
        for (const std::size_t bit : expected)
        {
            if (bit / 64 == wordIndex)
            {
                expectedWord |= 1ULL << (bit % 64);
            }
        }
        assert(selected.activeWord(wordIndex) == expectedWord);
    }

    std::vector<std::size_t> actual;
    for (
        std::size_t bit = selected.findFirst();
        bit < selected.size();
        bit = selected.findNext(bit)
    )
    {
        actual.push_back(bit);
    }
    assert(actual == expected);

    Mask copied = selected;
    assert(copied == selected);
    assert(std::hash<Mask>{}(copied) == std::hash<Mask>{}(selected));

    std::unordered_set<Mask> set{selected, copied};
    assert(set.size() == 1);
    assert(set.contains(selected));

    std::unordered_map<Mask, int> map;
    map.emplace(selected, 7);
    assert(map.at(copied) == 7);

    Mask odd;
    Mask even;
    for (std::size_t bit = 0; bit < width; bit++)
    {
        (bit % 2 == 0 ? even : odd).set(bit);
    }
    assert((odd & even).none());
    assert(odd.disjoint(even));
    assert(!odd.intersects(even));
    assert((odd | even) == full);
    assert(full.contains(odd));
    assert((full ^ odd) == even);
    assert(full.intersectionCount(odd) == odd.count());
    assert(selected.lowWordBelow(0) == 0);

    Mask compound = full;
    compound &= even;
    assert(compound == even);
    compound |= odd;
    assert(compound == full);
    compound ^= odd;
    assert(compound == even);

    if (width <= 64)
    {
        unsigned long long expectedWord = 0;
        for (const std::size_t bit : expected)
        {
            expectedWord |= 1ULL << bit;
        }
        assert(selected.to_ullong() == expectedWord);
        assert(selected.lowWordBelow(64) == expectedWord);
    }
    else
    {
        Mask high;
        high.set(64);
        bool overflowed = false;
        try
        {
            static_cast<void>(high.to_ullong());
        }
        catch (const std::overflow_error &)
        {
            overflowed = true;
        }
        assert(overflowed);
    }
}

int main()
{
    constexpr std::array<std::size_t, 25> widths{
        0, 1, 63, 64, 65, 127, 128, 129, 191, 192, 193, 255, 256,
        257, 319, 320, 321, 383, 384, 385, 447, 448, 449, 511, 512
    };
    for (const std::size_t width : widths)
    {
        testWidth<EdgeMask>(width);
        testWidth<AtomMask>(width);
    }

    EdgeMask::configure(65);
    AtomMask::configure(257);
    assert(EdgeMask::activeWordCount() == 2);
    assert(AtomMask::activeWordCount() == 5);

    bool rejected = false;
    try
    {
        EdgeMask::configure(513);
    }
    catch (const std::length_error &)
    {
        rejected = true;
    }
    assert(rejected);
    assert(EdgeMask::size() == 65);
}
