#pragma once

#include <cstddef>
#include <string_view>

/**
 * @brief Minimal UTF-8 decoding shared by the input parsers and the JSON writer.
 *
 * Both sides must agree on what a well-formed label is: the parsers reject text
 * they cannot represent, and the pathway writer escapes only text the parsers
 * already accepted. Keeping one decoder here stops the two from drifting.
 */
namespace utf8
{
    /** One decoded scalar value, or a failure when @ref length is zero. */
    struct DecodedCharacter
    {
        char32_t codePoint = 0;
        std::size_t length = 0;
    };

    /**
     * @brief Decode the sequence beginning at @p offset.
     *
     * Overlong encodings, surrogate halves, truncated sequences and values
     * above U+10FFFF all fail, so an accepted string survives a round trip
     * through any conforming JSON reader.
     *
     * @param text Bytes to decode
     * @param offset Index of the leading byte
     * @return The scalar value and its byte length, or a zero length on failure
     */
    [[nodiscard]] inline DecodedCharacter decode(
        std::string_view text,
        std::size_t offset
    )
    {
        if (offset >= text.size()) return {};

        const unsigned char leadingByte =
            static_cast<unsigned char>(text[offset]);
        if (leadingByte < 0x80) return {leadingByte, 1};

        std::size_t length = 0;
        char32_t codePoint = 0;
        if ((leadingByte & 0xE0) == 0xC0)
        {
            length = 2;
            codePoint = static_cast<char32_t>(leadingByte & 0x1F);
        }
        else if ((leadingByte & 0xF0) == 0xE0)
        {
            length = 3;
            codePoint = static_cast<char32_t>(leadingByte & 0x0F);
        }
        else if ((leadingByte & 0xF8) == 0xF0)
        {
            length = 4;
            codePoint = static_cast<char32_t>(leadingByte & 0x07);
        }
        else return {};

        if (offset + length > text.size()) return {};
        for (std::size_t index = 1; index < length; index++)
        {
            const unsigned char continuationByte =
                static_cast<unsigned char>(text[offset + index]);
            if ((continuationByte & 0xC0) != 0x80) return {};
            codePoint = (codePoint << 6) |
                static_cast<char32_t>(continuationByte & 0x3F);
        }

        static constexpr char32_t shortestForm[] = {0, 0, 0x80, 0x800, 0x10000};
        if (codePoint < shortestForm[length]) return {};
        if (codePoint > 0x10FFFF) return {};
        if (codePoint >= 0xD800 && codePoint <= 0xDFFF) return {};
        return {codePoint, length};
    }

    /**
     * @brief Report whether every byte of @p text belongs to a valid sequence.
     */
    [[nodiscard]] inline bool wellFormed(std::string_view text)
    {
        std::size_t offset = 0;
        while (offset < text.size())
        {
            const DecodedCharacter character = decode(text, offset);
            if (character.length == 0) return false;
            offset += character.length;
        }
        return true;
    }
}
