#pragma once

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

/**
 * @brief Runtime-width bit mask with an exact one-word small specialization.
 *
 * Domains up to 64 bits store their value directly in the eight-byte object.
 * Wider domains store an intrusive copy-on-write pointer in the same slot and
 * allocate exactly the configured number of words. Thus the common path has
 * the footprint and data movement of one uint64_t, while larger graphs have no
 * fixed mask ceiling and retain cheap value copies.
 *
 * A domain must be configured only after all masks from its preceding width
 * have been destroyed (persistent masks must be destroyed and reconstructed)
 * and all mask-keyed containers have been cleared. The surrounding search is
 * process-global and already observes that lifecycle; it is not re-entrant, so
 * reference counts are non-atomic.
 */
template<typename Domain>
class ActiveWordMask
{
public:
    using word_type = std::uint64_t;
    static constexpr std::size_t wordBits =
        std::numeric_limits<word_type>::digits;

    ActiveWordMask() noexcept
    {
        if (isSmall()) [[likely]] storage_.word = 0;
        else storage_.tail = nullptr;
    }

    ActiveWordMask(unsigned long long value)
    {
        if (isSmall()) [[likely]]
        {
            storage_.word = static_cast<word_type>(value) & firstWordMask();
        }
        else
        {
            storage_.tail = nullptr;
            if (value != 0)
            {
                storage_.tail = allocateWords();
                storage_.tail->data()[0] = static_cast<word_type>(value);
            }
        }
    }

    ActiveWordMask(const ActiveWordMask &other) noexcept
    {
        if (isSmall()) [[likely]] storage_.word = other.storage_.word;
        else
        {
            storage_.tail = other.storage_.tail;
            retain(storage_.tail);
        }
    }

    ActiveWordMask(ActiveWordMask &&other) noexcept
    {
        if (isSmall()) [[likely]]
        {
            storage_.word = std::exchange(other.storage_.word, word_type{0});
        }
        else
        {
            storage_.tail = std::exchange(other.storage_.tail, nullptr);
        }
    }

    [[gnu::always_inline]] ~ActiveWordMask()
    {
        if (!isSmall() && storage_.tail != nullptr) [[unlikely]]
        {
            destroyWide(storage_.tail);
        }
    }

    ActiveWordMask &operator=(const ActiveWordMask &other) noexcept
    {
        if (this == &other) return *this;
        if (isSmall()) [[likely]]
        {
            storage_.word = other.storage_.word;
            return *this;
        }
        retain(other.storage_.tail);
        release(storage_.tail);
        storage_.tail = other.storage_.tail;
        return *this;
    }

    ActiveWordMask &operator=(ActiveWordMask &&other) noexcept
    {
        if (this == &other) return *this;
        if (isSmall()) [[likely]]
        {
            storage_.word = std::exchange(other.storage_.word, word_type{0});
        }
        else
        {
            release(storage_.tail);
            storage_.tail = std::exchange(other.storage_.tail, nullptr);
        }
        return *this;
    }

    ActiveWordMask &operator=(unsigned long long value)
    {
        if (isSmall()) [[likely]]
        {
            storage_.word = static_cast<word_type>(value) & firstWordMask();
            return *this;
        }
        WideWords *replacement = nullptr;
        if (value != 0)
        {
            replacement = allocateWords();
            replacement->data()[0] = static_cast<word_type>(value);
        }
        release(storage_.tail);
        storage_.tail = replacement;
        return *this;
    }

    /** Configure this domain's logical width without a fixed bit ceiling. */
    static void configure(std::size_t bitCount) noexcept
    {
        clearArena();
        activeBitCount_ = bitCount;
        activeWordCount_ = bitCount / wordBits + (bitCount % wordBits != 0);
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
        return std::numeric_limits<std::size_t>::max();
    }

    [[nodiscard]] word_type activeWord(std::size_t index) const noexcept
    {
        if (isSmall()) [[likely]] return storage_.word;
        return storage_.tail == nullptr ? 0 : storage_.tail->data()[index];
    }

    [[gnu::always_inline]] ActiveWordMask &set(
        std::size_t position,
        bool value = true
    )
    {
        checkPosition(position);
        if (isSmall()) [[likely]]
        {
            const word_type bit = word_type{1} << position;
            if (value) storage_.word |= bit;
            else storage_.word &= ~bit;
            return *this;
        }
        const std::size_t wordIndex = position / wordBits;
        const word_type bit = word_type{1} << (position % wordBits);
        const word_type oldWord = activeWord(wordIndex);
        if (((oldWord & bit) != 0) == value) return *this;
        return setWideChanged(wordIndex, bit, value);
    }

    /** Return this mask with one bit set, avoiding a COW retain/detach pair. */
    [[nodiscard, gnu::always_inline]] ActiveWordMask withBitSet(
        std::size_t position
    ) const
    {
        checkPosition(position);
        if (isSmall()) [[likely]]
        {
            ActiveWordMask result;
            result.storage_.word = storage_.word | (word_type{1} << position);
            return result;
        }
        return withBitSetWide(position);
    }

    [[nodiscard, gnu::noinline]] ActiveWordMask withBitSetWide(
        std::size_t position
    ) const
    {
        const std::size_t wordIndex = position / wordBits;
        const word_type bit = word_type{1} << (position % wordBits);
        if ((activeWord(wordIndex) & bit) != 0) return *this;
        ActiveWordMask result;
        if (storage_.tail != nullptr)
        {
            result.storage_.tail = allocateWords(false);
            std::copy_n(
                storage_.tail->data(),
                activeWordCount_,
                result.storage_.tail->data()
            );
        }
        else result.storage_.tail = allocateWords();
        result.storage_.tail->data()[wordIndex] |= bit;
        return result;
    }

    /** Build a mask from exactly activeWordCount() logical words. */
    [[nodiscard]] static ActiveWordMask fromActiveWords(
        const word_type *words
    )
    {
        ActiveWordMask result;
        if (activeWordCount_ == 0) return result;
        if (isSmall()) [[likely]]
        {
            result.storage_.word = words[0] & lastWordMask();
            return result;
        }
        result.storage_.tail = allocateWords(false);
        std::copy_n(words, activeWordCount_, result.storage_.tail->data());
        result.storage_.tail->data()[activeWordCount_ - 1] &= lastWordMask();
        result.releaseWordsIfZero();
        return result;
    }

    ActiveWordMask &set()
    {
        if (activeWordCount_ <= 1) [[likely]]
        {
            storage_.word = lastWordMask();
            return *this;
        }
        WideWords *words = ensureUniqueWords();
        std::fill_n(
            words->data(),
            activeWordCount_,
            std::numeric_limits<word_type>::max()
        );
        words->data()[activeWordCount_ - 1] &= lastWordMask();
        return *this;
    }

    ActiveWordMask &reset() noexcept
    {
        if (isSmall()) [[likely]] storage_.word = 0;
        else
        {
            release(storage_.tail);
            storage_.tail = nullptr;
        }
        return *this;
    }

    ActiveWordMask &reset(std::size_t position)
    {
        return set(position, false);
    }

    ActiveWordMask &flip(std::size_t position)
    {
        checkPosition(position);
        if (isSmall()) [[likely]]
        {
            storage_.word ^= word_type{1} << position;
            return *this;
        }
        WideWords *words = ensureUniqueWords();
        words->data()[position / wordBits] ^=
            word_type{1} << (position % wordBits);
        releaseWordsIfZero();
        return *this;
    }

    ActiveWordMask &flip()
    {
        if (activeWordCount_ <= 1) [[likely]]
        {
            storage_.word = ~storage_.word & lastWordMask();
            return *this;
        }
        WideWords *words = ensureUniqueWords();
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            words->data()[i] = ~words->data()[i];
        }
        words->data()[activeWordCount_ - 1] &= lastWordMask();
        releaseWordsIfZero();
        return *this;
    }

    [[nodiscard]] bool operator[](std::size_t position) const noexcept
    {
        return (
            activeWord(position / wordBits) &
            (word_type{1} << (position % wordBits))
        ) != 0;
    }

    [[nodiscard]] bool test(std::size_t position) const
    {
        checkPosition(position);
        return (*this)[position];
    }

    [[nodiscard, gnu::always_inline]] std::size_t count() const noexcept
    {
        if (isSmall()) [[likely]]
        {
            return static_cast<std::size_t>(std::popcount(storage_.word));
        }
        return countWide();
    }

    [[nodiscard, gnu::noinline]] std::size_t countWide() const noexcept
    {
        if (storage_.tail == nullptr) return 0;
        std::size_t result = 0;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result += static_cast<std::size_t>(
                std::popcount(storage_.tail->data()[i])
            );
        }
        return result;
    }

    [[nodiscard]] bool any() const noexcept
    {
        if (isSmall()) [[likely]] return storage_.word != 0;
        return storage_.tail != nullptr;
    }

    [[nodiscard]] bool none() const noexcept
    {
        return !any();
    }

    [[nodiscard]] bool all() const noexcept
    {
        if (activeBitCount_ == 0) return true;
        if (isSmall()) [[likely]] return storage_.word == lastWordMask();
        if (storage_.tail == nullptr) return false;
        for (std::size_t i = 0; i + 1 < activeWordCount_; i++)
        {
            if (storage_.tail->data()[i] !=
                std::numeric_limits<word_type>::max()) return false;
        }
        return storage_.tail->data()[activeWordCount_ - 1] == lastWordMask();
    }

    [[nodiscard]] unsigned long long to_ullong() const
    {
        if (isSmall()) [[likely]] return storage_.word;
        for (std::size_t i = 1; i < activeWordCount_; i++)
        {
            if (activeWord(i) != 0)
            {
                throw std::overflow_error(
                    "mask does not fit in unsigned long long"
                );
            }
        }
        return activeWord(0);
    }

    [[nodiscard]] unsigned long long lowWordBelow(
        std::size_t limit
    ) const noexcept
    {
        if (activeWordCount_ == 0 || limit == 0) return 0;
        word_type word = activeWord(0);
        if (limit < wordBits) word &= (word_type{1} << limit) - 1;
        return word;
    }

    [[nodiscard, gnu::always_inline]] bool intersects(
        const ActiveWordMask &other
    ) const noexcept
    {
        if (isSmall()) [[likely]]
        {
            return (storage_.word & other.storage_.word) != 0;
        }
        return intersectsWide(other);
    }

    [[nodiscard, gnu::noinline]] bool intersectsWide(
        const ActiveWordMask &other
    ) const noexcept
    {
        if (storage_.tail == nullptr || other.storage_.tail == nullptr)
        {
            return false;
        }
        if (storage_.tail == other.storage_.tail) return true;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if ((storage_.tail->data()[i] & other.storage_.tail->data()[i]) != 0)
            {
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] bool disjoint(const ActiveWordMask &other) const noexcept
    {
        return !intersects(other);
    }

    [[nodiscard]] bool contains(const ActiveWordMask &other) const noexcept
    {
        if (isSmall()) [[likely]]
        {
            return (storage_.word & other.storage_.word) == other.storage_.word;
        }
        if (other.storage_.tail == nullptr) return true;
        if (storage_.tail == nullptr) return false;
        if (storage_.tail == other.storage_.tail) return true;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if ((storage_.tail->data()[i] & other.storage_.tail->data()[i]) !=
                other.storage_.tail->data()[i]) return false;
        }
        return true;
    }

    [[nodiscard]] std::size_t intersectionCount(
        const ActiveWordMask &other
    ) const noexcept
    {
        if (isSmall()) [[likely]]
        {
            return static_cast<std::size_t>(
                std::popcount(storage_.word & other.storage_.word)
            );
        }
        if (storage_.tail == nullptr || other.storage_.tail == nullptr) return 0;
        if (storage_.tail == other.storage_.tail) return count();
        std::size_t result = 0;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result += static_cast<std::size_t>(std::popcount(
                storage_.tail->data()[i] & other.storage_.tail->data()[i]
            ));
        }
        return result;
    }

    [[nodiscard]] std::size_t findFirst() const noexcept
    {
        if (isSmall()) [[likely]]
        {
            return storage_.word == 0
                ? activeBitCount_
                : static_cast<std::size_t>(std::countr_zero(storage_.word));
        }
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            const word_type word = activeWord(i);
            if (word != 0)
            {
                return i * wordBits +
                    static_cast<std::size_t>(std::countr_zero(word));
            }
        }
        return activeBitCount_;
    }

    [[nodiscard]] std::size_t findNext(std::size_t position) const noexcept
    {
        if (activeBitCount_ == 0 || position >= activeBitCount_ - 1)
        {
            return activeBitCount_;
        }
        position++;
        std::size_t wordIndex = position / wordBits;
        word_type word = activeWord(wordIndex);
        word &= std::numeric_limits<word_type>::max() << (position % wordBits);
        if (word != 0)
        {
            return wordIndex * wordBits +
                static_cast<std::size_t>(std::countr_zero(word));
        }
        for (wordIndex++; wordIndex < activeWordCount_; wordIndex++)
        {
            word = activeWord(wordIndex);
            if (word != 0)
            {
                return wordIndex * wordBits +
                    static_cast<std::size_t>(std::countr_zero(word));
            }
        }
        return activeBitCount_;
    }

    ActiveWordMask &operator&=(const ActiveWordMask &other)
    {
        if (isSmall()) [[likely]]
        {
            storage_.word &= other.storage_.word;
            return *this;
        }
        if (storage_.tail == nullptr) return *this;
        if (other.storage_.tail == nullptr)
        {
            return reset();
        }
        if (storage_.tail == other.storage_.tail) return *this;
        replaceWithComputed(storage_.tail, other.storage_.tail, Operation::and_);
        return *this;
    }

    [[gnu::always_inline]] ActiveWordMask &operator|=(
        const ActiveWordMask &other
    )
    {
        if (isSmall()) [[likely]]
        {
            storage_.word |= other.storage_.word;
            return *this;
        }
        return orAssignWide(other);
    }

    [[gnu::noinline]] ActiveWordMask &orAssignWide(
        const ActiveWordMask &other
    )
    {
        if (other.storage_.tail == nullptr ||
            storage_.tail == other.storage_.tail) return *this;
        if (storage_.tail == nullptr)
        {
            storage_.tail = other.storage_.tail;
            retain(storage_.tail);
            return *this;
        }
        replaceWithOr(storage_.tail, other.storage_.tail);
        return *this;
    }

    [[gnu::always_inline]] ActiveWordMask &operator^=(
        const ActiveWordMask &other
    )
    {
        if (isSmall()) [[likely]]
        {
            storage_.word ^= other.storage_.word;
            return *this;
        }
        return xorAssignWide(other);
    }

    [[gnu::noinline]] ActiveWordMask &xorAssignWide(
        const ActiveWordMask &other
    )
    {
        if (other.storage_.tail == nullptr) return *this;
        if (storage_.tail == nullptr)
        {
            storage_.tail = other.storage_.tail;
            retain(storage_.tail);
            return *this;
        }
        if (storage_.tail == other.storage_.tail) return reset();
        replaceWithComputed(storage_.tail, other.storage_.tail, Operation::xor_);
        return *this;
    }

    friend ActiveWordMask operator&(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    )
    {
        ActiveWordMask result;
        if (isSmall()) [[likely]]
        {
            result.storage_.word = left.storage_.word & right.storage_.word;
            return result;
        }
        if (left.storage_.tail == nullptr || right.storage_.tail == nullptr)
        {
            return result;
        }
        if (left.storage_.tail == right.storage_.tail)
        {
            result.storage_.tail = left.storage_.tail;
            retain(result.storage_.tail);
            return result;
        }
        result.storage_.tail = makeComputed(
            left.storage_.tail, right.storage_.tail, Operation::and_
        );
        return result;
    }

    friend ActiveWordMask operator|(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    )
    {
        ActiveWordMask result;
        if (isSmall()) [[likely]]
        {
            result.storage_.word = left.storage_.word | right.storage_.word;
            return result;
        }
        if (left.storage_.tail == nullptr ||
            left.storage_.tail == right.storage_.tail)
        {
            result.storage_.tail = right.storage_.tail;
            retain(result.storage_.tail);
            return result;
        }
        if (right.storage_.tail == nullptr)
        {
            result.storage_.tail = left.storage_.tail;
            retain(result.storage_.tail);
            return result;
        }
        result.storage_.tail = makeOr(left.storage_.tail, right.storage_.tail);
        return result;
    }

    friend ActiveWordMask operator^(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    )
    {
        ActiveWordMask result;
        if (isSmall()) [[likely]]
        {
            result.storage_.word = left.storage_.word ^ right.storage_.word;
            return result;
        }
        if (left.storage_.tail == right.storage_.tail) return result;
        if (left.storage_.tail == nullptr || right.storage_.tail == nullptr)
        {
            result.storage_.tail = left.storage_.tail == nullptr
                ? right.storage_.tail : left.storage_.tail;
            retain(result.storage_.tail);
            return result;
        }
        result.storage_.tail = makeComputed(
            left.storage_.tail, right.storage_.tail, Operation::xor_
        );
        return result;
    }

    friend ActiveWordMask operator~(const ActiveWordMask &mask)
    {
        ActiveWordMask result;
        if (activeWordCount_ == 0) return result;
        if (isSmall()) [[likely]]
        {
            result.storage_.word = ~mask.storage_.word & lastWordMask();
            return result;
        }
        result.storage_.tail = allocateWords(false);
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result.storage_.tail->data()[i] = ~mask.activeWord(i);
        }
        result.storage_.tail->data()[activeWordCount_ - 1] &= lastWordMask();
        result.releaseWordsIfZero();
        return result;
    }

    [[gnu::always_inline]] friend bool operator==(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    ) noexcept
    {
        if (isSmall()) [[likely]] return left.storage_.word == right.storage_.word;
        return equalsWide(left, right);
    }

    [[gnu::noinline]] static bool equalsWide(
        const ActiveWordMask &left,
        const ActiveWordMask &right
    ) noexcept
    {
        if (left.storage_.tail == right.storage_.tail) return true;
        if (left.storage_.tail == nullptr || right.storage_.tail == nullptr)
            return false;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if (left.storage_.tail->data()[i] != right.storage_.tail->data()[i])
            {
                return false;
            }
        }
        return true;
    }

    [[gnu::always_inline]] friend bool operator==(
        const ActiveWordMask &mask,
        unsigned long long value
    ) noexcept
    {
        const word_type trimmed = static_cast<word_type>(value) & firstWordMask();
        if (isSmall()) [[likely]] return mask.storage_.word == trimmed;
        if (trimmed == 0) return mask.storage_.tail == nullptr;
        if (mask.storage_.tail == nullptr) return false;
        if (mask.activeWord(0) != trimmed) return false;
        for (std::size_t i = 1; i < activeWordCount_; i++)
        {
            if (mask.activeWord(i) != 0) return false;
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

    [[nodiscard, gnu::always_inline]] std::size_t hash() const noexcept
    {
        if (activeWordCount_ == 0) return 0;
        if (isSmall()) [[likely]]
        {
            return std::hash<word_type>{}(storage_.word);
        }
        return hashWide();
    }

    [[nodiscard, gnu::noinline]] std::size_t hashWide() const noexcept
    {
        std::size_t result = std::hash<word_type>{}(activeWord(0));
        for (std::size_t i = 1; i < activeWordCount_; i++)
        {
            result ^= std::hash<word_type>{}(activeWord(i)) +
                static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) +
                (result << 6) + (result >> 2);
        }
        return result;
    }

private:
    struct WideWords
    {
        union
        {
            std::size_t references;
            WideWords *nextFree;
        };

        [[nodiscard]] word_type *data() noexcept
        {
            return reinterpret_cast<word_type *>(this + 1);
        }

        [[nodiscard]] const word_type *data() const noexcept
        {
            return reinterpret_cast<const word_type *>(this + 1);
        }
    };

    static_assert(sizeof(WideWords) == sizeof(std::size_t));

    struct alignas(std::max_align_t) ArenaBlock
    {
        ArenaBlock *next;
        std::size_t used;
        std::size_t capacity;

        [[nodiscard]] std::byte *data() noexcept
        {
            return reinterpret_cast<std::byte *>(this + 1);
        }
    };

    struct ArenaState
    {
        ArenaBlock *blocks = nullptr;
        WideWords *freeWords = nullptr;
    };

    union Storage
    {
        word_type word;
        WideWords *tail;

        constexpr Storage() noexcept: word(0) {}
    };

    enum class Operation { and_, or_, xor_ };

    [[nodiscard]] static bool isSmall() noexcept
    {
        return activeWordCount_ <= 1;
    }

    static void checkPosition(std::size_t position)
    {
        if (position >= activeBitCount_)
        {
            throw std::out_of_range("mask position outside active width");
        }
    }

    static constexpr word_type lowBits(std::size_t bitCount) noexcept
    {
        if (bitCount == 0) return 0;
        if (bitCount == wordBits)
        {
            return std::numeric_limits<word_type>::max();
        }
        return (word_type{1} << bitCount) - 1;
    }

    static word_type firstWordMask() noexcept
    {
        return lowBits(std::min(activeBitCount_, wordBits));
    }

    static word_type lastWordMask() noexcept
    {
        if (activeWordCount_ == 0) return 0;
        const std::size_t remainder = activeBitCount_ % wordBits;
        return lowBits(remainder == 0 ? wordBits : remainder);
    }

    static WideWords *allocateWords(bool initialise = true)
    {
        constexpr std::size_t headerSize = sizeof(WideWords);
        if (activeWordCount_ >
            (std::numeric_limits<std::size_t>::max() - headerSize) /
                sizeof(word_type))
        {
            throw std::bad_array_new_length();
        }
        const std::size_t rawSize =
            headerSize + activeWordCount_ * sizeof(word_type);
        constexpr std::size_t alignment = alignof(WideWords);
        if (rawSize > std::numeric_limits<std::size_t>::max() - (alignment - 1))
        {
            throw std::bad_array_new_length();
        }
        const std::size_t allocationSize = alignUp(rawSize);
        ArenaState &state = arena_;
        void *memory = nullptr;
        if (state.freeWords != nullptr)
        {
            WideWords *words = state.freeWords;
            state.freeWords = words->nextFree;
            memory = words;
        }
        else
        {
            ArenaBlock *block = state.blocks;
            if (block == nullptr ||
                allocationSize > block->capacity - block->used)
            {
                constexpr std::size_t defaultBlockSize = 64 * 1024;
                const std::size_t capacity = std::max(
                    defaultBlockSize, allocationSize
                );
                constexpr std::size_t blockHeader = sizeof(ArenaBlock);
                if (capacity >
                    std::numeric_limits<std::size_t>::max() - blockHeader)
                {
                    throw std::bad_array_new_length();
                }
                block = static_cast<ArenaBlock *>(
                    ::operator new(blockHeader + capacity)
                );
                block->next = state.blocks;
                block->used = 0;
                block->capacity = capacity;
                state.blocks = block;
            }
            memory = block->data() + block->used;
            block->used += allocationSize;
        }
        WideWords *words = static_cast<WideWords *>(memory);
        words->references = 1;
        if (initialise)
        {
            std::fill_n(words->data(), activeWordCount_, word_type{0});
        }
        return words;
    }

    static void retain(WideWords *words) noexcept
    {
        if (words == nullptr) return;
        if (words->references != std::numeric_limits<std::size_t>::max())
        {
            ++words->references;
        }
    }

    static void release(WideWords *words) noexcept
    {
        if (words == nullptr ||
            words->references == std::numeric_limits<std::size_t>::max()) return;
        if (--words->references == 0)
        {
            words->nextFree = arena_.freeWords;
            arena_.freeWords = words;
        }
    }

    [[gnu::noinline]] static void destroyWide(WideWords *words) noexcept
    {
        release(words);
    }

    static constexpr std::size_t alignUp(std::size_t value) noexcept
    {
        constexpr std::size_t alignment = alignof(WideWords);
        return (value + alignment - 1) & ~(alignment - 1);
    }

    static void clearArena() noexcept
    {
        ArenaBlock *block = arena_.blocks;
        while (block != nullptr)
        {
            ArenaBlock *next = block->next;
            ::operator delete(block);
            block = next;
        }
        arena_.blocks = nullptr;
        arena_.freeWords = nullptr;
    }

    WideWords *ensureUniqueWords()
    {
        if (storage_.tail == nullptr)
        {
            storage_.tail = allocateWords();
            return storage_.tail;
        }
        if (storage_.tail->references == 1)
        {
            return storage_.tail;
        }
        WideWords *replacement = allocateWords(false);
        std::copy_n(
            storage_.tail->data(),
            activeWordCount_,
            replacement->data()
        );
        release(storage_.tail);
        storage_.tail = replacement;
        return storage_.tail;
    }

    [[gnu::noinline]] ActiveWordMask &setWideChanged(
        std::size_t wordIndex,
        word_type bit,
        bool value
    )
    {
        WideWords *words = ensureUniqueWords();
        if (value) words->data()[wordIndex] |= bit;
        else words->data()[wordIndex] &= ~bit;
        if (!value) releaseWordsIfZero();
        return *this;
    }

    static bool wordsAreZero(const WideWords *words) noexcept
    {
        if (words == nullptr) return true;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            if (words->data()[i] != 0) return false;
        }
        return true;
    }

    void releaseWordsIfZero() noexcept
    {
        if (!wordsAreZero(storage_.tail)) return;
        release(storage_.tail);
        storage_.tail = nullptr;
    }

    static word_type apply(
        word_type left,
        word_type right,
        Operation operation
    ) noexcept
    {
        switch (operation)
        {
            case Operation::and_: return left & right;
            case Operation::or_: return left | right;
            case Operation::xor_: return left ^ right;
        }
        return 0;
    }

    static WideWords *makeComputed(
        const WideWords *left,
        const WideWords *right,
        Operation operation
    )
    {
        bool any = false;
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            any |= apply(left->data()[i], right->data()[i], operation) != 0;
        }
        if (!any) return nullptr;
        WideWords *result = allocateWords(false);
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result->data()[i] = apply(
                left->data()[i], right->data()[i], operation
            );
        }
        return result;
    }

    static WideWords *makeOr(
        const WideWords *left,
        const WideWords *right
    )
    {
        WideWords *result = allocateWords(false);
        for (std::size_t i = 0; i < activeWordCount_; i++)
        {
            result->data()[i] = left->data()[i] | right->data()[i];
        }
        return result;
    }

    void replaceWithOr(
        const WideWords *left,
        const WideWords *right
    )
    {
        if (storage_.tail->references == 1)
        {
            for (std::size_t i = 0; i < activeWordCount_; i++)
            {
                storage_.tail->data()[i] =
                    left->data()[i] | right->data()[i];
            }
            return;
        }
        WideWords *replacement = makeOr(left, right);
        release(storage_.tail);
        storage_.tail = replacement;
    }

    void replaceWithComputed(
        const WideWords *left,
        const WideWords *right,
        Operation operation
    )
    {
        if (storage_.tail->references == 1)
        {
            for (std::size_t i = 0; i < activeWordCount_; i++)
            {
                storage_.tail->data()[i] = apply(
                    left->data()[i], right->data()[i], operation
                );
            }
            releaseWordsIfZero();
            return;
        }
        WideWords *replacement = makeComputed(left, right, operation);
        release(storage_.tail);
        storage_.tail = replacement;
    }

    inline static std::size_t activeBitCount_ = wordBits;
    inline static std::size_t activeWordCount_ = 1;
    inline static ArenaState arena_;
    Storage storage_;
};

struct EdgeMaskDomain;
struct AtomMaskDomain;

using EdgeMask = ActiveWordMask<EdgeMaskDomain>;
using AtomMask = ActiveWordMask<AtomMaskDomain>;

static_assert(!std::is_same_v<EdgeMask, AtomMask>);
static_assert(sizeof(EdgeMask) == sizeof(std::uint64_t));
static_assert(sizeof(AtomMask) == sizeof(std::uint64_t));
static_assert(std::is_nothrow_copy_constructible_v<EdgeMask>);
static_assert(std::is_nothrow_move_constructible_v<EdgeMask>);
static_assert(std::is_nothrow_move_assignable_v<EdgeMask>);

namespace std
{
    template<typename Domain>
    struct hash<ActiveWordMask<Domain>>
    {
        size_t operator()(const ActiveWordMask<Domain> &mask) const noexcept
        {
            return mask.hash();
        }
    };
}
