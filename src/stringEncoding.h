#pragma once

#include <cstddef>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "utf8.h"

namespace parallelassemblycpp::detail::stringEncoding
{

/** Decode UTF-8 symbols, optionally appending their byte offsets and the end. */
inline std::u32string decodeInput(
    std::string_view input,
    std::vector<size_t> *byteOffsets = nullptr
)
{
    std::u32string result;
    for (size_t offset = 0; offset < input.size();)
    {
        const utf8::DecodedCharacter decoded = utf8::decode(input, offset);
        if (decoded.length == 0)
            throw std::invalid_argument(
                "string input is not valid UTF-8 at byte " + std::to_string(offset)
            );
        if (result.size() >= static_cast<size_t>(std::numeric_limits<int>::max()) - 1)
            throw std::invalid_argument("string is too long to index");
        if (byteOffsets != nullptr) byteOffsets->push_back(offset);
        result.push_back(decoded.codePoint);
        offset += decoded.length;
    }
    if (byteOffsets != nullptr) byteOffsets->push_back(input.size());
    return result;
}

/**
 * @brief Write one JSON string, escaping control bytes and all non-ASCII text.
 *
 * Non-ASCII text is decoded and re-emitted as \\uXXXX escapes rather than raw
 * bytes, so the pathway file is pure ASCII and stays readable whatever encoding
 * the consumer opens it with. Text that is not valid UTF-8 has no JSON
 * spelling at all, so it is reported instead of written.
 *
 * @throws std::runtime_error when @p value is not valid UTF-8
 */
inline void writeJsonString(std::string_view value, std::ostream &output)
{
    static constexpr char hexDigits[] = "0123456789ABCDEF";
    const auto writeUnitEscape = [&output](char32_t unit)
    {
        output << "\\u"
               << hexDigits[(unit >> 12) & 0x0f]
               << hexDigits[(unit >> 8) & 0x0f]
               << hexDigits[(unit >> 4) & 0x0f]
               << hexDigits[unit & 0x0f];
    };

    output.put('"');
    std::size_t offset = 0;
    while (offset < value.size())
    {
        const unsigned char character =
            static_cast<unsigned char>(value[offset]);
        if (character < 0x80)
        {
            offset++;
            switch (character)
            {
                case '"': output << "\\\""; break;
                case '\\': output << "\\\\"; break;
                case '\b': output << "\\b"; break;
                case '\f': output << "\\f"; break;
                case '\n': output << "\\n"; break;
                case '\r': output << "\\r"; break;
                case '\t': output << "\\t"; break;
                default:
                    if (character < 0x20) writeUnitEscape(character);
                    else output.put(static_cast<char>(character));
                    break;
            }
            continue;
        }

        const utf8::DecodedCharacter decoded = utf8::decode(value, offset);
        if (decoded.length == 0)
        {
            throw std::runtime_error(
                "cannot write JSON: text is not valid UTF-8"
            );
        }
        offset += decoded.length;

        if (decoded.codePoint < 0x10000) writeUnitEscape(decoded.codePoint);
        else
        {
            const char32_t remainder = decoded.codePoint - 0x10000;
            writeUnitEscape(0xD800 + (remainder >> 10));
            writeUnitEscape(0xDC00 + (remainder & 0x3FF));
        }
    }
    output.put('"');
}

} // namespace parallelassemblycpp::detail::stringEncoding
