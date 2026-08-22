// Compile this file directly to exercise the pathway JSON helpers without
// writing output files or running an assembly search.
#define ASSEMBLYCPP_NO_MAIN
#include "../v5/main.cpp"

#include <climits>
#include <cstdlib>
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

void testBondColoursAreAlwaysJsonValues()
{
    requireEqual(bondColour(SHRT_MIN), "\"error\"");
    requireEqual(bondColour(-1), "\"error\"");
    requireEqual(bondColour(0), "\"error\"");
    requireEqual(bondColour(1), "\"single\"");
    requireEqual(bondColour(2), "\"double\"");
    requireEqual(bondColour(3), "\"triple\"");
    requireEqual(bondColour(4), "\"4\"");
    requireEqual(
        bondColour(SHRT_MAX),
        "\"" + to_string(static_cast<int>(SHRT_MAX)) + "\""
    );
}

int main()
{
    testJsonStringEscaping();
    testBondColoursAreAlwaysJsonValues();
    return 0;
}
