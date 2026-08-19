#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>

inline constexpr std::size_t MASK_BIT_CAPACITY = 512;

/**
 * @brief Fixed-capacity bit mask whose work is proportional to its active words.
 *
 * Each mask domain is configured once per calculation before masks are inserted
 * into hash tables. Storage stays inline, while construction, copying, bitwise
 * operations, equality, population counts, set-bit searches, and hashing touch
 * only the words required by the current molecule.
 */
template<typename Domain, std::size_t Capacity = MASK_BIT_CAPACITY>
class ActiveWordMask
{
public:
    using word_type = std::uint64_t;
    static constexpr std::size_t wordBits =
        std::numeric_limits<word_type>::digits;
    static constexpr std::size_t wordCapacity =
        (Capacity + wordBits - 1) / wordBits;

    ActiveWordMask() noexcept
    {
        clearActiveWords();
    }

    ActiveWordMask(unsigned long long value) noexcept
    {
        assign(value);
    }

    ActiveWordMask(const ActiveWordMask &other) noexcept
    {
        copyActiveWords(other);
    }

    ActiveWordMask(ActiveWordMask &&other) noexcept
    {
        copyActiveWords(other);
    }

    ActiveWordMask &operator=(const ActiveWordMask &other) noexcept
    {
        if (this != &other) copyActiveWords(other);
        return *this;
    }

    ActiveWordMask &operator=(ActiveWordMask &&other) noexcept
    {
        if (this != &other) copyActiveWords(other);
        return *this;
    }

    ActiveWordMask &operator=(unsigned long long value) noexcept
    {
        assign(value);
        return *this;
    }

    /**
     * @brief Set the logical width for this mask domain.
     *
     * Reconfiguration invalidates the meaning of existing masks and therefore
     * must happen only after domain-specific hash containers have been cleared.
     */
    static void configure(std::size_t bitCount)
    {
        if (bitCount > Capacity)
        {
            throw std::length_error("mask capacity exceeded");
        }
        activeBitCount_ = bitCount;
        activeWordCount_ = (bitCount + wordBits - 1) / wordBits;
    }

    [[nodiscard]] static std::size_t size() noexcept
    {
        return activeBitCount_;
    }

    [[nodiscard]] static std::size_t activeWordCount() noexcept
    {
        return activeWordCount_;
    }

    [[nodiscard]] static constexpr std::size_t capacity() noexcept
    {
        return Capacity;
    }

    /**
     * @brief Read one word inside the configured mask width.
     *
     * The caller must pass an index below activeWordCount(). This deliberately
     * does not expose inactive inline storage, whose value is unspecified.
     */
    [[nodiscard]] word_type activeWord(std::size_t index) const noexcept
    {
        return words_[index];
    }

    ActiveWordMask &set(std::size_t position, bool value = true)
    {
        checkPosition(position);
        const word_type bit = word_type{1} << (position % wordBits);
        word_type &word = words_[position / wordBits];
        if (value) word |= bit;
        else word &= ~bit;
        return *this;
    }

    ActiveWordMask &set() noexcept
    {
        if (activeWordCount_ == 1)
        {
            words_[0] = lastWordMask();
            return *this;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            words_[i] = std::numeric_limits<word_type>::max();
        }
        trimLastWord();
        return *this;
    }

    ActiveWordMask &reset() noexcept
    {
        clearActiveWords();
        return *this;
    }

    ActiveWordMask &reset(std::size_t position)
    {
        return set(position, false);
    }

    ActiveWordMask &flip(std::size_t position)
    {
        checkPosition(position);
        words_[position / wordBits] ^=
            word_type{1} << (position % wordBits);
        return *this;
    }

    ActiveWordMask &flip() noexcept
    {
        if (activeWordCount_ == 1)
        {
            words_[0] = ~words_[0] & lastWordMask();
            return *this;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++) words_[i] = ~words_[i];
        trimLastWord();
        return *this;
    }

    [[nodiscard]] bool operator[](std::size_t position) const noexcept
    {
        return (
            words_[position / wordBits] &
            (word_type{1} << (position % wordBits))
        ) != 0;
    }

    [[nodiscard]] bool test(std::size_t position) const
    {
        checkPosition(position);
        return (*this)[position];
    }

    [[nodiscard]] std::size_t count() const noexcept
    {
        if (activeWordCount_ == 1)
        {
            return static_cast<std::size_t>(std::popcount(words_[0]));
        }
        std::size_t result = 0;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result += static_cast<std::size_t>(std::popcount(words_[i]));
        }
        return result;
    }

    [[nodiscard]] bool any() const noexcept
    {
        if (activeWordCount_ == 1) return words_[0] != 0;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if (words_[i] != 0) return true;
        }
        return false;
    }

    [[nodiscard]] bool none() const noexcept
    {
        return !any();
    }

    [[nodiscard]] bool all() const noexcept
    {
        if (activeBitCount_ == 0) return true;
        const std::size_t last = activeWordCount_ - 1;
        for (std::size_t i = 0; i < last; i++)
        {
            if (words_[i] != std::numeric_limits<word_type>::max()) return false;
        }
        return words_[last] == lastWordMask();
    }

    [[nodiscard]] unsigned long long to_ullong() const
    {
        for (std::size_t i = 1; i < activeWordCount_; i++)
        {
            if (words_[i] != 0)
            {
                throw std::overflow_error("mask does not fit in unsigned long long");
            }
        }
        return activeWordCount_ == 0 ? 0 : words_[0];
    }

    [[nodiscard]] unsigned long long lowWordBelow(
        std::size_t limit
    ) const noexcept
    {
        if (activeWordCount_ == 0 || limit == 0) return 0;
        word_type word = words_[0];
        if (limit < wordBits) word &= (word_type{1} << limit) - 1;
        return word;
    }

    [[nodiscard]] bool intersects(const ActiveWordMask &other) const noexcept
    {
        if (activeWordCount_ == 1)
        {
            return (words_[0] & other.words_[0]) != 0;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if ((words_[i] & other.words_[i]) != 0) return true;
        }
        return false;
    }

    [[nodiscard]] bool disjoint(const ActiveWordMask &other) const noexcept
    {
        return !intersects(other);
    }

    [[nodiscard]] bool contains(const ActiveWordMask &other) const noexcept
    {
        if (activeWordCount_ == 1)
        {
            return (words_[0] & other.words_[0]) == other.words_[0];
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if ((words_[i] & other.words_[i]) != other.words_[i]) return false;
        }
        return true;
    }

    [[nodiscard]] std::size_t intersectionCount(
        const ActiveWordMask &other
    ) const noexcept
    {
        if (activeWordCount_ == 1)
        {
            return static_cast<std::size_t>(
                std::popcount(words_[0] & other.words_[0])
            );
        }
        std::size_t result = 0;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result += static_cast<std::size_t>(
                std::popcount(words_[i] & other.words_[i])
            );
        }
        return result;
    }

    [[nodiscard]] std::size_t findFirst() const noexcept
    {
        if (activeWordCount_ == 1)
        {
            return words_[0] == 0
                ? activeBitCount_
                : static_cast<std::size_t>(std::countr_zero(words_[0]));
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if (words_[i] != 0)
            {
                return i * wordBits +
                    static_cast<std::size_t>(std::countr_zero(words_[i]));
            }
        }
        return activeBitCount_;
    }

    [[nodiscard]] std::size_t findNext(std::size_t position) const noexcept
    {
        if (position >= activeBitCount_ - (activeBitCount_ != 0))
        {
            return activeBitCount_;
        }

        position++;
        std::size_t wordIndex = position / wordBits;
        word_type word = words_[wordIndex];
        word &= std::numeric_limits<word_type>::max() << (position % wordBits);
        if (word != 0)
        {
            return wordIndex * wordBits +
                static_cast<std::size_t>(std::countr_zero(word));
        }

        for (wordIndex++; wordIndex < activeWordCount_; wordIndex++)
        {
            if (words_[wordIndex] != 0)
            {
                return wordIndex * wordBits + static_cast<std::size_t>(
                    std::countr_zero(words_[wordIndex])
                );
            }
        }
        return activeBitCount_;
    }

    ActiveWordMask &operator&=(const ActiveWordMask &other) noexcept
    {
        if (activeWordCount_ == 1)
        {
            words_[0] &= other.words_[0];
            return *this;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            words_[i] &= other.words_[i];
        }
        return *this;
    }

    ActiveWordMask &operator|=(const ActiveWordMask &other) noexcept
    {
        if (activeWordCount_ == 1)
        {
            words_[0] |= other.words_[0];
            return *this;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            words_[i] |= other.words_[i];
        }
        return *this;
    }

    ActiveWordMask &operator^=(const ActiveWordMask &other) noexcept
    {
        if (activeWordCount_ == 1)
        {
            words_[0] ^= other.words_[0];
            return *this;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            words_[i] ^= other.words_[i];
        }
        return *this;
    }

    friend ActiveWordMask operator&(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    ) noexcept
    {
        ActiveWordMask result(uninitialised);
        if (activeWordCount_ == 1)
        {
            result.words_[0] = left.words_[0] & right.words_[0];
            return result;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result.words_[i] = left.words_[i] & right.words_[i];
        }
        return result;
    }

    friend ActiveWordMask operator|(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    ) noexcept
    {
        ActiveWordMask result(uninitialised);
        if (activeWordCount_ == 1)
        {
            result.words_[0] = left.words_[0] | right.words_[0];
            return result;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result.words_[i] = left.words_[i] | right.words_[i];
        }
        return result;
    }

    friend ActiveWordMask operator^(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    ) noexcept
    {
        ActiveWordMask result(uninitialised);
        if (activeWordCount_ == 1)
        {
            result.words_[0] = left.words_[0] ^ right.words_[0];
            return result;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result.words_[i] = left.words_[i] ^ right.words_[i];
        }
        return result;
    }

    friend ActiveWordMask operator~(const ActiveWordMask &mask) noexcept
    {
        ActiveWordMask result(uninitialised);
        if (activeWordCount_ == 1)
        {
            result.words_[0] = ~mask.words_[0] & lastWordMask();
            return result;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result.words_[i] = ~mask.words_[i];
        }
        result.trimLastWord();
        return result;
    }

    friend bool operator==(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    ) noexcept
    {
        if (activeWordCount_ == 1) return left.words_[0] == right.words_[0];
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if (left.words_[i] != right.words_[i]) return false;
        }
        return true;
    }

    friend bool operator==(
        const ActiveWordMask &mask,
        unsigned long long value
    ) noexcept
    {
        if (activeWordCount_ == 0) return true;
        const word_type trimmed = static_cast<word_type>(value) & firstWordMask();
        if (mask.words_[0] != trimmed) return false;
        for (std::size_t i = 1; i < activeWordCount_; i++)
        {
            if (mask.words_[i] != 0) return false;
        }
        return true;
    }

    friend bool operator==(
        unsigned long long value,
        const ActiveWordMask &mask
    ) noexcept
    {
        return mask == value;
    }

    [[nodiscard]] std::size_t hash() const noexcept
    {
        if (activeWordCount_ == 0) return 0;
        if (activeWordCount_ == 1)
        {
            return std::hash<word_type>{}(words_[0]);
        }

        std::size_t result = std::hash<word_type>{}(words_[0]);
        for (std::size_t i = 1; i < activeWordCount_; i++)
        {
            result ^= std::hash<word_type>{}(words_[i]) +
                static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) +
                (result << 6) + (result >> 2);
        }
        return result;
    }

private:
    struct Uninitialised {};
    inline static constexpr Uninitialised uninitialised{};

    explicit ActiveWordMask(Uninitialised) noexcept {}

    static void checkPosition(std::size_t position)
    {
        if (position >= activeBitCount_)
        {
            throw std::out_of_range("mask position outside active width");
        }
    }

    static constexpr word_type lowBits(std::size_t bitCount) noexcept
    {
        return bitCount == 0 || bitCount == wordBits
            ? std::numeric_limits<word_type>::max()
            : (word_type{1} << bitCount) - 1;
    }

    static word_type firstWordMask() noexcept
    {
        return lowBits(activeBitCount_ < wordBits ? activeBitCount_ : wordBits);
    }

    static word_type lastWordMask() noexcept
    {
        const std::size_t remainder = activeBitCount_ % wordBits;
        return lowBits(remainder == 0 ? wordBits : remainder);
    }

    void trimLastWord() noexcept
    {
        if (activeWordCount_ != 0)
        {
            words_[activeWordCount_ - 1] &= lastWordMask();
        }
    }

    void clearActiveWords() noexcept
    {
        if (activeWordCount_ == 1)
        {
            words_[0] = 0;
            return;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++) words_[i] = 0;
    }

    void copyActiveWords(const ActiveWordMask &other) noexcept
    {
        if (activeWordCount_ == 1)
        {
            words_[0] = other.words_[0];
            return;
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            words_[i] = other.words_[i];
        }
    }

    void assign(unsigned long long value) noexcept
    {
        if (activeWordCount_ == 0) return;
        words_[0] = static_cast<word_type>(value) & firstWordMask();
        for (std::size_t i = 1; i < activeWordCount_; i++) words_[i] = 0;
    }

    inline static std::size_t activeBitCount_ = Capacity;
    inline static std::size_t activeWordCount_ = wordCapacity;
    word_type words_[wordCapacity];
};

struct EdgeMaskDomain;
struct AtomMaskDomain;

using EdgeMask = ActiveWordMask<EdgeMaskDomain>;
using AtomMask = ActiveWordMask<AtomMaskDomain>;

static_assert(!std::is_same_v<EdgeMask, AtomMask>);
static_assert(sizeof(EdgeMask) == MASK_BIT_CAPACITY / 8);
static_assert(sizeof(AtomMask) == MASK_BIT_CAPACITY / 8);

namespace std
{
    template<typename Domain, size_t Capacity>
    struct hash<ActiveWordMask<Domain, Capacity>>
    {
        size_t operator()(
            const ActiveWordMask<Domain, Capacity> &mask
        ) const noexcept
        {
            return mask.hash();
        }
    };
}
