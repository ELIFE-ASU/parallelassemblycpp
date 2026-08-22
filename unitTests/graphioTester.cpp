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
            left.mg.size() != right.mg.size()
        ) return false;

        for (size_t vertex = 0; vertex < left.mg.size(); vertex++)
        {
            if (
                left.mg[vertex].type != right.mg[vertex].type ||
                left.mg[vertex].list.size() != right.mg[vertex].list.size()
            ) return false;
            for (size_t edge = 0; edge < left.mg[vertex].list.size(); edge++)
            {
                const bond &leftBond = left.mg[vertex].list[edge];
                const bond &rightBond = right.mg[vertex].list[edge];
                if (
                    leftBond.n != rightBond.n ||
                    leftBond.type != rightBond.type
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
    verbose = false;

    molGraph graph = parse(validGraph);
    assert(graph.mg.size() == 3);
    assert(graph.totalBonds == 2);
    assert(graph.mg[0].type == "C");
    assert(graph.mg[1].type == "N");
    assert(graph.mg[2].type == "O");
    assert(graph.mg[0].list.size() == 1);
    assert(graph.mg[1].list.size() == 2);
    assert(graph.mg[2].list.size() == 1);
    assert(graph.mg[0].list[0].n == 1);
    assert(graph.mg[1].list[0].type == 1);
    assert(graph.mg[1].list[1].type == 2);

    molGraph replacementTarget = sentinelGraph();
    istringstream validInput(validGraph);
    graphio(validInput, replacementTarget);
    assert(sameGraph(replacementTarget, graph));

    const molGraph empty = parse("empty\n0\n\n\n\n");
    assert(empty.mg.empty());
    assert(empty.totalBonds == 0);

    const molGraph explicitlyPositive = parse(
        "positive signs\n+2\n+1 +2\nC C\n+1\n"
    );
    assert(explicitlyPositive.mg.size() == 2);
    assert(explicitlyPositive.totalBonds == 1);
    assert(explicitlyPositive.mg[0].list[0].type == 1);

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
        "exactly one integer"
    );
    expectRejected(
        "bad count\n-1\n\n\n\n",
        "short-index range"
    );
    expectRejected(
        "bad count\n" +
            to_string(static_cast<long long>(numeric_limits<short>::max()) + 1) +
            "\n\n\n\n",
        "short-index range"
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
        "expected 1 bond labels, found 2"
    );
    expectRejected(
        "negative type\n2\n1 2\nC C\n-1\n",
        "nonnegative short range"
    );
    expectRejected(
        "wide type\n2\n1 2\nC C\n" +
            to_string(static_cast<long long>(numeric_limits<short>::max()) + 1) +
            "\n",
        "nonnegative short range"
    );
    expectRejected(
        "invalid type\n2\n1 2\nC C\n1.5\n",
        "bond label must be an integer"
    );

    molGraph longestSupportedPath;
    const size_t vertexCount = static_cast<size_t>(numeric_limits<short>::max());
    longestSupportedPath.mg.reserve(vertexCount);
    for (size_t vertex = 0; vertex < vertexCount; vertex++)
        longestSupportedPath.addAtom("C");
    for (size_t vertex = 1; vertex < vertexCount; vertex++)
        longestSupportedPath.addBond(vertex - 1, vertex, 1);
    assert(longestSupportedPath.disjointFragments() == 1);

    return 0;
}
