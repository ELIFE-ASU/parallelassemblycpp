// Compile this file directly to exercise pathway witness ownership and JSON
// helpers without writing output files or running an assembly search.
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

void testRetainedPathwaySurvivesDecisionOwners(size_t edgeCount)
{
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(edgeCount);
    std::construct_at(std::addressof(allEdges));

    assemblyPathWitness witness;
    {
        // Initial root decisions borrow masks owned by enumeration; deeper
        // decisions borrow raw immutable words retained in the runtime DAG.
        EdgeMask rootMatch;
        EdgeMask rootDuplicate;
        rootMatch.set(0);
        rootDuplicate.set(edgeCount - 1);
        vector<uint64_t> dagMatch(EdgeMask::activeWordCount(), 0);
        vector<uint64_t> dagDuplicate(EdgeMask::activeWordCount(), 0);
        dagMatch.at(1 / EdgeMask::wordBits) |=
            uint64_t{1} << (1 % EdgeMask::wordBits);
        dagDuplicate.at((edgeCount - 2) / EdgeMask::wordBits) |=
            uint64_t{1} << ((edgeCount - 2) % EdgeMask::wordBits);

        witness.pushDecision(rootMatch, rootDuplicate);
        witness.pushDecision(
            EdgeMaskView::fromWords(dagMatch.data()),
            EdgeMaskView::fromWords(dagDuplicate.data())
        );
        witness.retainCurrent();
        witness.clearDecisions();

        // Reuse both kinds of producer storage before destroying them. The
        // winning witness must own its words rather than retaining views.
        rootMatch.reset();
        rootDuplicate.reset();
        fill(dagMatch.begin(), dagMatch.end(), 0);
        fill(dagDuplicate.begin(), dagDuplicate.end(), 0);
    }
    if (
        witness.best.size() != 2 ||
        witness.best[0].match.count() != 1 ||
        !witness.best[0].match[0] ||
        witness.best[0].duplicate.count() != 1 ||
        !witness.best[0].duplicate[edgeCount - 1] ||
        witness.best[1].match.count() != 1 ||
        !witness.best[1].match[1] ||
        witness.best[1].duplicate.count() != 1 ||
        !witness.best[1].duplicate[edgeCount - 2]
    ) abort();

    // A later improvement at a shorter depth replaces the whole checkpoint.
    EdgeMask replacement;
    replacement.set(edgeCount / 2);
    witness.pushDecision(replacement, replacement);
    witness.retainCurrent();
    witness.clearDecisions();
    if (
        witness.best.size() != 1 ||
        witness.best.front().match != replacement ||
        witness.best.front().duplicate != replacement
    ) abort();

    // The empty root witness must clear a previously retained pathway.
    witness.retainCurrent();
    if (!witness.best.empty()) abort();
}

void testPathwayCheckpointDivergentBranches()
{
    std::destroy_at(std::addressof(allEdges));
    EdgeMask::configure(129);
    std::construct_at(std::addressof(allEdges));

    assemblyPathWitness witness;
    EdgeMask first;
    EdgeMask second;
    EdgeMask third;
    first.set(0);
    second.set(64);
    third.set(128);

    witness.pushDecision(first, second);
    witness.retainCurrent();
    witness.pushDecision(second, third);
    witness.retainCurrent();
    witness.pushDecision(third, first);
    witness.retainCurrent();
    if (witness.best.size() != 3) abort();

    // Returning to depth one invalidates both saved descendants. A sibling
    // must retain the common ancestor and replace the entire old suffix.
    witness.popDecision();
    witness.popDecision();
    witness.pushDecision(third, second);
    witness.retainCurrent();
    if (
        witness.best.size() != 2 ||
        witness.best[0].match != first ||
        witness.best[0].duplicate != second ||
        witness.best[1].match != third ||
        witness.best[1].duplicate != second
    ) abort();

    // A shorter winning branch removes the previous suffix without changing
    // its surviving prefix; a new root then invalidates that prefix too.
    witness.popDecision();
    witness.retainCurrent();
    if (witness.best.size() != 1 || witness.best[0].match != first) abort();
    witness.popDecision();
    witness.pushDecision(second, first);
    witness.retainCurrent();
    witness.popDecision();
    if (
        witness.best.size() != 1 ||
        witness.best[0].match != second ||
        witness.best[0].duplicate != first
    ) abort();
}

int main()
{
    testJsonStringEscaping();
    testNonAsciiJsonStringEscaping();
    testBondColoursAreAlwaysJsonValues();
    testRetainedPathwaySurvivesDecisionOwners(32);
    testRetainedPathwaySurvivesDecisionOwners(129);
    testPathwayCheckpointDivergentBranches();
    return 0;
}
