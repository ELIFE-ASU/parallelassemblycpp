// This standalone test keeps its checks active in release-style builds.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <csignal>
#include <ctime>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
using IntegerVector = vector<int>;
using BooleanVector = vector<bool>;
using IntegerPair = pair<int, int>;

#include "../v5/activeWordMask.h"

constexpr int ceilLog2(int value)
{
    return std::bit_width(static_cast<unsigned int>(value - 1));
}

#include "../v5/globalPrimitives.h"
#include "../v5/ufds.h"
#include "../v5/molGraph.h"
#include "../v5/graphio.h"

namespace
{
    const string validGraph =
        "native graph with a descriptive name\n"
        "+3\n"
        "  1 2   2 3  \n"
        "C N O\n"
        "+1 2\n";

    molGraph sentinelGraph()
    {
        molGraph graph;
        string carbon = "sentinel-C";
        string oxygen = "sentinel-O";
        graph.addAtom(carbon);
        graph.addAtom(oxygen);
        graph.addBond(0, 1, 7);
        return graph;
    }

    bool sameGraph(const molGraph &left, const molGraph &right)
    {
        if (
            left.totalBonds != right.totalBonds ||
            left.atoms.size() != right.atoms.size()
        ) return false;

        for (size_t vertex = 0; vertex < left.atoms.size(); vertex++)
        {
            if (
                left.atoms[vertex].atomType !=
                    right.atoms[vertex].atomType ||
                left.atoms[vertex].bonds.size() !=
                    right.atoms[vertex].bonds.size()
            ) return false;
            for (
                size_t edge = 0;
                edge < left.atoms[vertex].bonds.size();
                edge++
            )
            {
                const bond &leftBond = left.atoms[vertex].bonds[edge];
                const bond &rightBond = right.atoms[vertex].bonds[edge];
                if (
                    leftBond.neighbourAtomIndex !=
                        rightBond.neighbourAtomIndex ||
                    leftBond.bondType != rightBond.bondType
                ) return false;
            }
        }
        return true;
    }

    molGraph parse(const string &source)
    {
        istringstream input(source);
        molGraph graph;
        graphio(input, graph);
        return graph;
    }

    void expectRejected(const string &source, const string &diagnosticFragment)
    {
        molGraph destination = sentinelGraph();
        const molGraph original = destination;
        istringstream input(source);
        bool rejected = false;
        try
        {
            graphio(input, destination);
        }
        catch (const runtime_error &error)
        {
            rejected = string(error.what()).find(diagnosticFragment) != string::npos;
        }
        assert(rejected);
        assert(sameGraph(destination, original));
    }
}

int main()
{
    removeHydrogens = true;
    verbose = false;

    molGraph graph = parse(validGraph);
    assert(graph.atoms.size() == 3);
    assert(graph.totalBonds == 2);
    assert(graph.atoms[0].atomType == "C");
    assert(graph.atoms[1].atomType == "N");
    assert(graph.atoms[2].atomType == "O");
    assert(graph.atoms[0].bonds.size() == 1);
    assert(graph.atoms[1].bonds.size() == 2);
    assert(graph.atoms[2].bonds.size() == 1);
    assert(graph.atoms[0].bonds[0].neighbourAtomIndex == 1);
    assert(graph.atoms[1].bonds[0].bondType == 1);
    assert(graph.atoms[1].bonds[1].bondType == 2);

    const string explicitHydrogenGraph =
        "explicit hydrogens\n"
        "4\n"
        "1 3 2 3 3 4\n"
        "H H C C\n"
        "1 1 1\n";
    const molGraph hydrogensRemoved = parse(explicitHydrogenGraph);
    assert(hydrogensRemoved.atoms.size() == 2);
    assert(hydrogensRemoved.totalBonds == 1);
    assert(hydrogensRemoved.atoms[0].atomType == "C");
    assert(hydrogensRemoved.atoms[1].atomType == "C");
    assert(hydrogensRemoved.atoms[0].bonds[0].neighbourAtomIndex == 1);

    removeHydrogens = false;
    const molGraph hydrogensRetained = parse(explicitHydrogenGraph);
    assert(hydrogensRetained.atoms.size() == 4);
    assert(hydrogensRetained.totalBonds == 3);
    assert(hydrogensRetained.atoms[0].atomType == "H");
    assert(hydrogensRetained.atoms[1].atomType == "H");

    removeHydrogens = true;
    const molGraph unrestrictedLabels = parse(
        "unrestricted labels\n"
        "2\n"
        "1 2\n"
        "COLLAPSE He\n"
        "1\n"
    );
    assert(unrestrictedLabels.atoms.size() == 2);
    assert(unrestrictedLabels.totalBonds == 1);
    assert(unrestrictedLabels.atoms[0].atomType == "COLLAPSE");
    assert(unrestrictedLabels.atoms[1].atomType == "He");

    molGraph replacementTarget = sentinelGraph();
    istringstream validInput(validGraph);
    graphio(validInput, replacementTarget);
    assert(sameGraph(replacementTarget, graph));

    const molGraph empty = parse("empty\n0\n\n\n\n");
    assert(empty.atoms.empty());
    assert(empty.totalBonds == 0);

    const molGraph explicitlyPositive = parse(
        "positive signs\n+2\n+1 +2\nC C\n+1\n"
    );
    assert(explicitlyPositive.atoms.size() == 2);
    assert(explicitlyPositive.totalBonds == 1);
    assert(explicitlyPositive.atoms[0].bonds[0].bondType == 1);

    expectRejected(
        "truncated\n2\n1 2\nC C\n",
        "missing bond label line"
    );
    expectRejected(
        "bad count\nnot-an-integer\n\n\n\n",
        "graph size must be an integer"
    );
    expectRejected(
        "bad count\n2 3\n\nC C\n\n",
        "must contain one integer"
    );
    expectRejected(
        "bad count\n-1\n\n\n\n",
        "graph size must be between 0 and"
    );
    expectRejected(
        "bad count\n" +
            to_string(static_cast<long long>(numeric_limits<short>::max()) + 1) +
            "\n\n\n\n",
        "graph size must be between 0 and"
    );
    expectRejected(
        "odd endpoints\n2\n1\nC C\n1\n",
        "complete endpoint pairs"
    );
    expectRejected(
        "out of range\n2\n1 3\nC C\n1\n",
        "outside the declared graph size"
    );
    expectRejected(
        "self loop\n2\n1 1\nC C\n1\n",
        "self-loop"
    );
    expectRejected(
        "missing atom\n2\n1 2\nC\n1\n",
        "expected 2 atom labels, found 1"
    );
    expectRejected(
        "extra atom\n2\n1 2\nC C O\n1\n",
        "expected 2 atom labels, found 3"
    );
    expectRejected(
        "missing type\n3\n1 2 2 3\nC C C\n1\n",
        "expected 2 bond labels, found 1"
    );
    expectRejected(
        "extra type\n2\n1 2\nC C\n1 2\n",
        "expected 1 bond label, found 2"
    );
    expectRejected(
        "zero type\n2\n1 2\nC C\n0\n",
        "zero-order bonds are not supported"
    );
    expectRejected(
        "negative type\n2\n1 2\nC C\n-1\n",
        "bond label must be between 0 and"
    );
    expectRejected(
        "wide type\n2\n1 2\nC C\n" +
            to_string(static_cast<long long>(numeric_limits<short>::max()) + 1) +
            "\n",
        "bond label must be between 0 and"
    );
    expectRejected(
        "invalid type\n2\n1 2\nC C\n1.5\n",
        "bond label must be an integer"
    );

    molGraph longestSupportedPath;
    const size_t vertexCount = static_cast<size_t>(numeric_limits<short>::max());
    longestSupportedPath.atoms.reserve(vertexCount);
    for (size_t vertex = 0; vertex < vertexCount; vertex++)
        longestSupportedPath.addAtom("C");
    for (size_t vertex = 1; vertex < vertexCount; vertex++)
    {
        const int secondVertex = static_cast<int>(vertex);
        longestSupportedPath.addBond(secondVertex - 1, secondVertex, 1);
    }
    assert(longestSupportedPath.disjointFragments() == 1);

    return 0;
}
