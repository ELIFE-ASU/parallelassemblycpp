// Compile this file directly to exercise the pathway JSON helpers without
// writing output files or running an assembly search.
#define PARALLELASSEMBLYCPP_NO_MAIN
#include "../src/main.cpp"

#include <cstdlib>
#include <limits>
#include <sstream>

void requireEqual(const string &actual, const string &expected)
{
    if (actual != expected) abort();
}

string jsonString(const string &value)
{
    ostringstream output;
    printJsonString(value, output);
    return output.str();
}

string bondColour(short type)
{
    ostringstream output;
    printBondColour(type, output);
    return output.str();
}

void testJsonStringEscaping()
{
    requireEqual(jsonString("C"), "\"C\"");
    requireEqual(jsonString("C\\\"N"), "\"C\\\\\\\"N\"");
    requireEqual(
        jsonString(string{"\b\f\n\r\t"}),
        "\"\\b\\f\\n\\r\\t\""
    );

    static constexpr char hexDigits[] = "0123456789ABCDEF";
    for (unsigned int value = 0; value < 0x20; value++)
    {
        string expected = "\"";
        switch (value)
        {
            case '\b': expected += "\\b"; break;
            case '\f': expected += "\\f"; break;
            case '\n': expected += "\\n"; break;
            case '\r': expected += "\\r"; break;
            case '\t': expected += "\\t"; break;
            default:
                expected += "\\u00";
                expected += hexDigits[value >> 4];
                expected += hexDigits[value & 0x0f];
                break;
        }
        expected += '"';
        requireEqual(jsonString(string(1, static_cast<char>(value))), expected);
    }
}

// Atom labels used to reach the pathway file as raw bytes, so a mis-encoded
// molfile produced JSON that no reader could decode.
void testNonAsciiJsonStringEscaping()
{
    requireEqual(jsonString("\xc3\x85"), "\"\\u00C5\"");
    requireEqual(jsonString("\xe2\x88\x9e"), "\"\\u221E\"");
    requireEqual(jsonString("\xf0\x9f\x98\x80"), "\"\\uD83D\\uDE00\"");
    requireEqual(jsonString("C\xc3\x85N"), "\"C\\u00C5N\"");

    const string malformed[] = {
        "\x80",
        "\xc3",
        "\xc3\x28",
        "\xc0\xaf",
        "\xed\xa0\x80",
        "\xf5\x80\x80\x80"
    };
    for (const string &value : malformed)
    {
        bool rejected = false;
        try
        {
            static_cast<void>(jsonString(value));
        }
        catch (const runtime_error &)
        {
            rejected = true;
        }
        if (!rejected) abort();
    }
}

void testBondColoursAreAlwaysJsonValues()
{
    requireEqual(
        bondColour(numeric_limits<short>::min()),
        "\"error\""
    );
    requireEqual(bondColour(-1), "\"error\"");
    requireEqual(bondColour(0), "\"error\"");
    requireEqual(bondColour(1), "\"single\"");
    requireEqual(bondColour(2), "\"double\"");
    requireEqual(bondColour(3), "\"triple\"");
    requireEqual(bondColour(4), "\"4\"");
    requireEqual(
        bondColour(numeric_limits<short>::max()),
        "\"" + to_string(numeric_limits<short>::max()) + "\""
    );
}

int main()
{
    testJsonStringEscaping();
    testNonAsciiJsonStringEscaping();
    testBondColoursAreAlwaysJsonValues();
    return 0;
}
