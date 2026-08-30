#include <algorithm>
#include <array>
#include <bit>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cstdint>
#include <csignal>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;
using vi = vector<int>;
using vb = vector<bool>;
using pii = pair<int, int>;

#include "../v5/activeWordMask.h"

constexpr int ceilLog2(int value)
{
    return std::bit_width(static_cast<unsigned int>(value - 1));
}

#include "../v5/globalPrimitives.h"
#include "../v5/ufds.h"
#include "../v5/molGraph.h"
#include "../v5/molfileParser.h"
#include "../v5/ioflag.h"

namespace
{
    const string validMolfile =
        "Fixed field test\n"
        "AssemblyCpp parser test\n"
        "\n"
        "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.0000    0.0000    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.0000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    4.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  2  0  0  0  0\n"
        "  3  4  3  0  0  0  0\n"
        "  4  5  1  0  0  0  0\n"
        "M  END\n";

    molGraph parse(const string &source)
    {
        istringstream input(source);
        molGraph graph;
        molfileParser(input, graph);
        return graph;
    }
}

int main(int argc, char **argv)
{
    removeHydrogens = false;
    verbose = false;

    if (argc > 1)
    {
        for (int argument = 1; argument < argc; argument++)
        {
            ifstream input(argv[argument]);
            if (!input)
            {
                cerr << "could not open " << argv[argument] << '\n';
                return 1;
            }
            molGraph graph;
            try
            {
                molfileParser(input, graph);
            }
            catch (const exception &error)
            {
                cerr << argv[argument] << ": " << error.what() << '\n';
                return 1;
            }
        }
        return 0;
    }

    ostringstream diagnostics;
    streambuf *const originalOutput = cout.rdbuf(diagnostics.rdbuf());
    molGraph graph = parse(validMolfile);
    cout.rdbuf(originalOutput);

    assert(diagnostics.str().empty());
    assert(graph.mg.size() == 5);
    assert(graph.totalBonds == 4);
    assert(graph.mg[0].type == "C");
    assert(graph.mg[1].type == "Cl");
    assert(graph.mg[2].type == "Si");
    assert(graph.mg[3].type == "Br");
    assert(graph.mg[4].type == "H");
    assert(graph.btypeS(1, 1) == 2);
    assert(graph.btypeS(2, 1) == 3);

    verbose = true;
    diagnostics.str("");
    diagnostics.clear();
    cout.rdbuf(diagnostics.rdbuf());
    static_cast<void>(parse(validMolfile));
    cout.rdbuf(originalOutput);
    assert(diagnostics.str().find("Molfile: 5 atoms, 4 bonds") != string::npos);
    assert(diagnostics.str().find("Graph: 5 atoms, 4 bonds") != string::npos);

    verbose = false;
    string invalidBond = validMolfile;
    const size_t bond = invalidBond.find("  4  5  1");
    assert(bond != string::npos);
    invalidBond.replace(bond, 9, "  4  6  1");
    bool rejected = false;
    try
    {
        static_cast<void>(parse(invalidBond));
    }
    catch (const runtime_error &error)
    {
        rejected = string(error.what()).find(
            "bond endpoint is outside the atom range"
        ) != string::npos;
    }
    assert(rejected);

    string selfBond = validMolfile;
    const size_t selfBondPosition = selfBond.find("  4  5  1");
    assert(selfBondPosition != string::npos);
    selfBond.replace(selfBondPosition, 9, "  4  4  1");
    rejected = false;
    try
    {
        static_cast<void>(parse(selfBond));
    }
    catch (const runtime_error &error)
    {
        rejected = string(error.what()).find("self-loop") != string::npos;
    }
    assert(rejected);

    string zeroOrderBond = validMolfile;
    const size_t zeroOrderPosition = zeroOrderBond.find("  4  5  1");
    assert(zeroOrderPosition != string::npos);
    zeroOrderBond.replace(zeroOrderPosition, 9, "  4  5  0");
    molGraph zeroOrderDestination;
    string zeroOrderSentinel = "zero-order-sentinel";
    zeroOrderDestination.addAtom(zeroOrderSentinel);
    istringstream zeroOrderInput(zeroOrderBond);
    rejected = false;
    try
    {
        molfileParser(zeroOrderInput, zeroOrderDestination);
    }
    catch (const runtime_error &error)
    {
        rejected = string(error.what()).find(
            "zero-order bonds are not supported"
        ) != string::npos;
    }
    assert(rejected);
    assert(zeroOrderDestination.mg.size() == 1);
    assert(zeroOrderDestination.mg.front().type == zeroOrderSentinel);
    assert(zeroOrderDestination.totalBonds == 0);

    string negativeBondOrder = validMolfile;
    const size_t negativeBondPosition = negativeBondOrder.find("  4  5  1");
    assert(negativeBondPosition != string::npos);
    negativeBondOrder.replace(negativeBondPosition, 9, "  4  5 -1");
    rejected = false;
    try
    {
        static_cast<void>(parse(negativeBondOrder));
    }
    catch (const runtime_error &error)
    {
        rejected = string(error.what()).find(
            "bond order must be nonnegative"
        ) != string::npos;
    }
    assert(rejected);

    molGraph unchanged;
    string sentinel = "sentinel";
    unchanged.addAtom(sentinel);
    istringstream invalidInput(invalidBond);
    try
    {
        molfileParser(invalidInput, unchanged);
    }
    catch (const runtime_error &)
    {
    }
    assert(unchanged.mg.size() == 1);
    assert(unchanged.mg.front().type == sentinel);

    string unsupportedVersion = validMolfile;
    const size_t version = unsupportedVersion.find("V2000");
    assert(version != string::npos);
    unsupportedVersion.replace(version, 5, "V3000");
    rejected = false;
    try
    {
        static_cast<void>(parse(unsupportedVersion));
    }
    catch (const runtime_error &error)
    {
        rejected = string(error.what()).find("expected a V2000") != string::npos;
    }
    assert(rejected);

    char executable[] = "AssemblyCpp";
    char option[] = "--verbose=1";
    char input[] = "input.mol";
    char *arguments[] = {executable, option, input};
    const CommandLineArguments parsed = parseCommandLine(3, arguments);
    assert(parsed.input == input);
    assert(verbose);

    return 0;
}
