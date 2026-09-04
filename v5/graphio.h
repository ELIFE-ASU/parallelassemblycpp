#pragma once

#include <charconv>
#include <iostream>
#include <istream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace graphioDetail
{
    inline void readLine(
        std::istream &input,
        std::string &line,
        const char *description
    )
    {
        if (!std::getline(input, line))
        {
            throw std::runtime_error(
                std::string("invalid native graph: missing ") + description
            );
        }
    }

    [[nodiscard]] inline std::vector<std::string> tokens(
        const std::string &line
    )
    {
        std::istringstream input(line);
        std::vector<std::string> result;
        for (std::string token; input >> token; )
            result.push_back(std::move(token));
        return result;
    }

    [[nodiscard]] inline long long integerToken(
        std::string_view token,
        const char *description
    )
    {
        // The legacy parser accepted an explicit positive sign, so preserve it
        // when using from_chars.
        if (!token.empty() && token.front() == '+') token.remove_prefix(1);
        long long value = 0;
        const char *begin = token.data();
        const char *end = begin + token.size();
        const auto conversion = std::from_chars(begin, end, value);
        if (
            token.empty() ||
            conversion.ec != std::errc() ||
            conversion.ptr != end
        )
        {
            throw std::runtime_error(
                std::string("invalid native graph: ") + description +
                " must be an integer"
            );
        }
        return value;
    }

    [[nodiscard]] inline std::size_t graphSize(const std::string &line)
    {
        const std::vector<std::string> fields = tokens(line);
        if (fields.size() != 1)
        {
            throw std::runtime_error(
                "invalid native graph: graph size line must contain one integer"
            );
        }

        const long long value = integerToken(fields.front(), "graph size");
        if (
            value < 0 ||
            value > static_cast<long long>(std::numeric_limits<short>::max())
        )
        {
            throw std::runtime_error(
                "invalid native graph: graph size must be between 0 and " +
                std::to_string(std::numeric_limits<short>::max())
            );
        }
        return static_cast<std::size_t>(value);
    }

    [[nodiscard]] inline std::vector<std::pair<int, int>> edges(
        const std::string &line,
        std::size_t vertexCount
    )
    {
        const std::vector<std::string> fields = tokens(line);
        if (fields.size() % 2 != 0)
        {
            throw std::runtime_error(
                "invalid native graph: edge list must contain "
                "complete endpoint pairs"
            );
        }

        std::vector<std::pair<int, int>> result;
        result.reserve(fields.size() / 2);
        std::vector<std::size_t> degrees(vertexCount, 0);
        constexpr std::size_t maximumDegree =
            static_cast<std::size_t>(std::numeric_limits<short>::max()) + 1;
        for (std::size_t field = 0; field < fields.size(); field += 2)
        {
            const long long first = integerToken(
                fields[field], "edge endpoint"
            );
            const long long second = integerToken(
                fields[field + 1], "edge endpoint"
            );
            if (
                first < 1 || second < 1 ||
                first > static_cast<long long>(vertexCount) ||
                second > static_cast<long long>(vertexCount)
            )
            {
                throw std::runtime_error(
                    "invalid native graph: edge endpoint is outside the "
                    "declared graph size"
                );
            }
            if (first == second)
            {
                throw std::runtime_error(
                    "invalid native graph: self-loop edges are not supported"
                );
            }
            const std::size_t firstIndex = static_cast<std::size_t>(first - 1);
            const std::size_t secondIndex = static_cast<std::size_t>(second - 1);
            if (
                ++degrees[firstIndex] > maximumDegree ||
                ++degrees[secondIndex] > maximumDegree
            )
            {
                throw std::runtime_error(
                    "invalid native graph: vertex degree exceeds the supported "
                    "maximum of " + std::to_string(maximumDegree)
                );
            }
            result.emplace_back(
                static_cast<int>(firstIndex),
                static_cast<int>(secondIndex)
            );
        }
        return result;
    }

    inline void requireCardinality(
        std::size_t actual,
        std::size_t expected,
        const char *singularDescription
    )
    {
        if (actual == expected) return;
        throw std::runtime_error(
            std::string("invalid native graph: expected ") +
            std::to_string(expected) + ' ' + singularDescription +
            (expected == 1 ? "" : "s") + ", found " +
            std::to_string(actual)
        );
    }
}

/**
 * @brief Parse AssemblyCpp's five-line native graph format transactionally.
 *
 * @param inputStream input stream
 * @param molecule output molGraph, replaced only after the complete input is valid
 */
inline void graphio(std::istream &inputStream, molGraph &molecule)
{
    std::string nameLine;
    std::string sizeLine;
    std::string endpointLine;
    std::string atomLine;
    std::string bondLine;
    graphioDetail::readLine(inputStream, nameLine, "graph name line");
    graphioDetail::readLine(inputStream, sizeLine, "graph size line");
    graphioDetail::readLine(inputStream, endpointLine, "edge endpoint line");
    graphioDetail::readLine(inputStream, atomLine, "atom label line");
    graphioDetail::readLine(inputStream, bondLine, "bond label line");

    const std::size_t vertexCount = graphioDetail::graphSize(sizeLine);
    const std::vector<std::pair<int, int>> edgeList =
        graphioDetail::edges(endpointLine, vertexCount);
    std::vector<std::string> atomLabels = graphioDetail::tokens(atomLine);
    const std::vector<std::string> bondFields = graphioDetail::tokens(bondLine);
    graphioDetail::requireCardinality(
        atomLabels.size(), vertexCount, "atom label"
    );
    graphioDetail::requireCardinality(
        bondFields.size(), edgeList.size(), "bond label"
    );

    std::vector<short> bondLabels;
    bondLabels.reserve(bondFields.size());
    for (const std::string &field : bondFields)
    {
        const long long value = graphioDetail::integerToken(field, "bond label");
        if (
            value < 0 ||
            value > static_cast<long long>(std::numeric_limits<short>::max())
        )
        {
            throw std::runtime_error(
                "invalid native graph: bond label must be between 0 and " +
                std::to_string(std::numeric_limits<short>::max())
            );
        }
        if (value == 0)
        {
            throw std::runtime_error(
                "invalid native graph: zero-order bonds are not supported"
            );
        }
        bondLabels.push_back(static_cast<short>(value));
    }

    molGraph parsed;
    parsed.atoms.reserve(vertexCount);
    for (std::string &label : atomLabels) parsed.addAtom(label);
    for (std::size_t edge = 0; edge < edgeList.size(); edge++)
    {
        parsed.addBond(
            edgeList[edge].first,
            edgeList[edge].second,
            bondLabels[edge]
        );
    }

    if (removeHydrogens) parsed.removeExplicitHydrogens();
    if (verbose)
    {
        const std::vector<std::string> nameFields = graphioDetail::tokens(nameLine);
        const std::string graphName = nameFields.empty() ? "" : nameFields.front();
        std::cout << "Native graph: " << graphName << '\n';
        parsed.printToCout();
    }
    molecule = std::move(parsed);
}
