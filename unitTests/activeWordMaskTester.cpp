// This executable is also built by Release/CI presets; keep its checks active.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../v5/activeWordMask.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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
static_assert(sizeof(EdgeMask) == sizeof(std::uint64_t));
static_assert(sizeof(AtomMask) == sizeof(std::uint64_t));

std::vector<std::size_t> boundaryBits(std::size_t width)
{
    std::vector<std::size_t> bits;
    if (width == 0) return bits;

    bits.push_back(0);
    bits.push_back(width - 1);
    for (std::size_t boundary = 64; boundary < width; boundary += 64)
    {
        bits.push_back(boundary - 1);
        bits.push_back(boundary);
        if (boundary + 1 < width) bits.push_back(boundary + 1);
    }
    std::ranges::sort(bits);
    const auto newEnd = std::ranges::unique(bits).begin();
    bits.erase(newEnd, bits.end());
    return bits;
}

template<typename Mask>
void testWidth(std::size_t width)
{
    Mask::configure(width);
    assert(Mask::size() == width);
    assert(Mask::activeWordCount() == (width + 63) / 64);
    assert(Mask::capacity() >= width);

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

    Mask selected;
    const std::vector<std::size_t> expected = boundaryBits(width);
    for (const std::size_t bit : expected)
    {
        selected.set(bit);
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

    // Parallel root jobs cross workers as plain words and are reconstructed
    // only after the receiving worker configures its own mask domain.
    std::vector<std::uint64_t> serialized(Mask::activeWordCount());
    for (std::size_t word = 0; word < serialized.size(); ++word)
        serialized[word] = selected.activeWord(word);
    Mask reconstructed = Mask::fromActiveWords(serialized.data());
    assert(reconstructed == selected);
    if (width != 0)
    {
        reconstructed.flip(width - 1);
        assert(reconstructed != selected);
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

    // A copied wide mask must detach before mutation. These checks apply to
    // small masks as well, preserving identical value semantics at 64/65.
    if (width != 0)
    {
        const std::size_t changedBit = width > 64 ? 64 : width - 1;
        const bool originalValue = selected.test(changedBit);
        copied.flip(changedBit);
        assert(copied.test(changedBit) != originalValue);
        assert(selected.test(changedBit) == originalValue);
        copied.flip(changedBit);
        assert(copied == selected);
    }

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

    bool rejectedBoundary = false;
    try
    {
        selected.set(width);
    }
    catch (const std::out_of_range &)
    {
        rejectedBoundary = true;
    }
    assert(rejectedBoundary);
}

template<typename Mask>
void testCopyMoveAndContainers(std::size_t width)
{
    Mask::configure(width);
    assert(width > 256);

    Mask original;
    for (const std::size_t bit : boundaryBits(width)) original.set(bit);

    Mask copyConstructed = original;
    Mask copyAssigned;
    copyAssigned = original;
    copyConstructed.flip(1);
    copyAssigned.reset(0);
    assert(original.test(1) != copyConstructed.test(1));
    assert(original.test(0));
    assert(!copyAssigned.test(0));

    Mask expected = original;
    Mask moveConstructed(std::move(copyConstructed));
    assert(moveConstructed.test(1) != expected.test(1));
    copyConstructed = expected;
    assert(copyConstructed == expected);

    Mask moveAssigned;
    moveAssigned.set(width - 2);
    moveAssigned = std::move(copyAssigned);
    assert(!moveAssigned.test(0));
    copyAssigned = expected;
    assert(copyAssigned == expected);

    // Reallocation repeatedly copies or moves masks which share their wide
    // storage. Mutating one element must not affect any other element.
    std::vector<Mask> copies;
    for (std::size_t index = 0; index < 257; index++)
        copies.push_back(original);
    assert(copies.front() == original);
    assert(copies.back() == original);
    copies[128].flip(2);
    assert(copies[128] != original);
    assert(copies[127] == original);
    assert(copies[129] == original);

    // Force several unordered-container rehashes with distinct wide keys.
    std::unordered_set<Mask> keys;
    keys.max_load_factor(0.25F);
    std::vector<Mask> expectedKeys;
    for (std::size_t bit = 3; bit < 131; bit++)
    {
        Mask key = original.withBitSet(bit);
        if (key == original) key.flip(bit);
        expectedKeys.push_back(key);
        keys.insert(key);
    }
    assert(keys.size() == expectedKeys.size());
    keys.rehash(2048);
    for (const Mask &key : expectedKeys) assert(keys.contains(key));

    std::unordered_map<Mask, std::size_t> values;
    values.max_load_factor(0.25F);
    for (std::size_t index = 0; index < expectedKeys.size(); index++)
    {
        values.emplace(expectedKeys[index], index);
    }
    values.rehash(2048);
    for (std::size_t index = 0; index < expectedKeys.size(); index++)
    {
        assert(values.at(expectedKeys[index]) == index);
    }
}

void testAccumulatorWidth(std::size_t width)
{
    EdgeMask::configure(width);
    const std::size_t wordCount = EdgeMask::activeWordCount();

    EdgeMask first;
    EdgeMask second;
    EdgeMask probe;
    EdgeMask endProbe;
    first.set(0).set(width - 1);
    second.set(width / 2).set(width - 1);
    probe.set(width - 1);
    endProbe.set(0).set(width - 1);
    if (width > 64) probe.set(64);
    const EdgeMask expected = first | second;

    EdgeMaskAccumulatorBuffer buffer(3);
    assert(buffer.size() == 3);
    assert(buffer.span().size() == 3);
    assert(buffer[0].activeWordCount() == wordCount);
    for (const EdgeMaskAccumulator &accumulator : buffer)
        assert(accumulator.none());

    buffer[0].add(first);
    buffer[1] |= second;
    buffer[2].add(buffer[0]).add(buffer[1]);
    assert(buffer[0].toEdgeMask() == first);
    assert(buffer[1].toEdgeMask() == second);
    assert(buffer[2].toEdgeMask() == expected);
    assert(buffer[2].any());
    assert(buffer[2].count() == expected.count());
    assert(
        buffer[2].intersectionCount(probe) ==
        expected.intersectionCount(probe)
    );
    assert(
        buffer[2].intersectionCount(endProbe) ==
        expected.intersectionCount(endProbe)
    );
    EdgeMask retained = expected;
    retained.intersectWords(buffer[0]);
    assert(retained == first);
    assert(expected == (first | second));
    for (std::size_t word = 0; word < wordCount; ++word)
        assert(buffer[2].activeWord(word) == expected.activeWord(word));

    const EdgeMaskAccumulatorBuffer &constBuffer = buffer;
    assert(constBuffer.span().size() == buffer.size());
    assert(constBuffer.span().data() == buffer.data());

    // Growing preserves the existing flat rows and clears only the new suffix.
    buffer.resize(7);
    assert(buffer[0].toEdgeMask() == first);
    assert(buffer[1].toEdgeMask() == second);
    assert(buffer[2].toEdgeMask() == expected);
    for (std::size_t index = 3; index < buffer.size(); ++index)
        assert(buffer[index].none());

    // Distinct wide accumulators must remain isolated after the flat word
    // vector and the view vector have both had an opportunity to reallocate.
    buffer[3] |= probe;
    buffer[4] |= first;
    assert(buffer[3].toEdgeMask() == probe);
    assert(buffer[4].toEdgeMask() == first);
    assert(buffer[0].toEdgeMask() == first);

    buffer.resetRange(1, 2);
    assert(buffer[0].toEdgeMask() == first);
    assert(buffer[1].none());
    assert(buffer[2].none());
    assert(buffer[3].toEdgeMask() == probe);

    // Shrinking and regrowing retains the prefix and zeroes re-exposed rows.
    buffer.resize(1);
    assert(buffer[0].toEdgeMask() == first);
    buffer.resize(5);
    assert(buffer[0].toEdgeMask() == first);
    for (std::size_t index = 1; index < buffer.size(); ++index)
        assert(buffer[index].none());

    buffer.reset();
    for (const EdgeMaskAccumulator &accumulator : buffer)
        assert(accumulator.none());

    buffer[0] |= expected;
    const std::size_t retainedCapacity = buffer.capacity();
    buffer.clear();
    assert(buffer.empty());
    assert(buffer.capacity() == retainedCapacity);
    buffer.reset(5);
    assert(buffer.capacity() == retainedCapacity);
    for (const EdgeMaskAccumulator &accumulator : buffer)
        assert(accumulator.none());

    bool rejectedRange = false;
    try
    {
        buffer.resetRange(buffer.size(), 1);
    }
    catch (const std::out_of_range &)
    {
        rejectedRange = true;
    }
    assert(rejectedRange);

    // Standalone accumulators are allocation-free for the inline routes.
    if (wordCount <= EdgeMaskAccumulator::inlineWordCapacity)
    {
        EdgeMaskAccumulator firstAccumulator;
        EdgeMaskAccumulator secondAccumulator;
        firstAccumulator.add(first);
        secondAccumulator.add(second);
        EdgeMaskAccumulator inlineAccumulator;
        inlineAccumulator |= firstAccumulator;
        inlineAccumulator |= secondAccumulator;
        assert(inlineAccumulator.toEdgeMask() == expected);
        assert(inlineAccumulator.intersectionCount(endProbe) == 2);
        inlineAccumulator.clear();
        assert(inlineAccumulator.none());
    }
}

void testDomainIndependence()
{
    EdgeMask::configure(1025);
    AtomMask::configure(65);
    {
        EdgeMask edges;
        AtomMask atoms;
        edges.set(1024).set(512).set(0);
        atoms.set(64).set(0);
        assert(edges.count() == 3);
        assert(atoms.count() == 2);
        assert(EdgeMask::size() == 1025);
        assert(AtomMask::size() == 65);
        assert(EdgeMask::activeWordCount() == 17);
        assert(AtomMask::activeWordCount() == 2);
    }

    // Configuring one domain cannot change the other domain's logical width.
    EdgeMask::configure(513);
    assert(EdgeMask::size() == 513);
    assert(AtomMask::size() == 65);
    assert(EdgeMask::activeWordCount() == 9);
    assert(AtomMask::activeWordCount() == 2);
}

template<typename Mask>
void testSequentialReconfiguration()
{
    constexpr std::array<std::size_t, 10> widths{
        0, 1, 64, 65, 512, 513, 1024, 1025, 63, 2049
    };
    for (const std::size_t width : widths)
    {
        Mask::configure(width);
        {
            Mask value;
            if (width != 0)
            {
                value.set(0).set(width - 1);
                assert(value.test(0));
                assert(value.test(width - 1));
                assert(value.count() == (width == 1 ? 1 : 2));
            }
            else
            {
                assert(value.none());
                assert(value.all());
            }
        }
    }

    // Return the domain to its allocation-free small representation after
    // validating a wide-to-small transition.
    Mask::configure(64);
    Mask value;
    value.set(63);
    assert(value.count() == 1);
}

int main()
{
    constexpr std::array<std::size_t, 33> widths{
        0, 1, 63, 64, 65, 127, 128, 129, 191, 192, 193, 255, 256,
        257, 319, 320, 321, 383, 384, 385, 447, 448, 449, 511, 512,
        513, 575, 576, 577, 1023, 1024, 1025, 2049
    };
    for (const std::size_t width : widths)
    {
        testWidth<EdgeMask>(width);
        testWidth<AtomMask>(width);
    }

    testCopyMoveAndContainers<EdgeMask>(1025);
    testCopyMoveAndContainers<AtomMask>(1025);
    constexpr std::array<std::size_t, 6> accumulatorWidths{
        64, 65, 127, 128, 129, 257
    };
    for (const std::size_t width : accumulatorWidths)
        testAccumulatorWidth(width);
    testDomainIndependence();
    testSequentialReconfiguration<EdgeMask>();
    testSequentialReconfiguration<AtomMask>();
}
