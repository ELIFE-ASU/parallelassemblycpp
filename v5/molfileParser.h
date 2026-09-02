#pragma once

#include <charconv>
#include <iostream>
#include <istream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace molfileParserDetail
{
    [[nodiscard]] inline std::string_view fixedField(
        const std::string &line,
        std::size_t offset,
        std::size_t width,
        const char *description
    )
    {
        if (line.size() < offset + width)
        {
            throw std::runtime_error(
                std::string("invalid molfile: ") + description + " is too short"
            );
        }
        return std::string_view(line).substr(offset, width);
    }

    [[nodiscard]] inline std::string_view trimSpaces(std::string_view field)
    {
        while (!field.empty() && field.front() == ' ') field.remove_prefix(1);
        while (!field.empty() && field.back() == ' ') field.remove_suffix(1);
        return field;
    }

    [[nodiscard]] inline int integerField(
        const std::string &line,
        std::size_t offset,
        std::size_t width,
        const char *description
    )
    {
        const std::string_view field = trimSpaces(
            fixedField(line, offset, width, description)
        );
        int value = 0;
        const auto result = std::from_chars(
            field.data(),
            field.data() + field.size(),
            value
        );
        if (
            field.empty() ||
            result.ec != std::errc() ||
            result.ptr != field.data() + field.size()
        )
        {
            throw std::runtime_error(
                std::string("invalid molfile: ") + description +
                " is not a valid integer"
            );
        }
        return value;
    }

    inline void readLine(
        std::istream &molfile,
        std::string &line,
        const char *description
    )
    {
        if (!std::getline(molfile, line))
        {
            throw std::runtime_error(
                std::string("invalid molfile: missing ") + description
            );
        }
    }
}

/** Parse a V2000 molfile into a molecular graph. */
void molfileParser(std::istream &molfile, molGraph &molecule)
{
    std::string currLine;
    molGraph parsed;

    for (int headerLine = 0; headerLine < 3; headerLine++)
    {
        molfileParserDetail::readLine(molfile, currLine, "header line");
    }

    molfileParserDetail::readLine(molfile, currLine, "counts line");
    if (
        currLine.size() < 39 ||
        std::string_view(currLine).substr(34, 5) != "V2000"
    )
    {
        throw std::runtime_error(
            "unsupported molfile format: expected a V2000 counts line"
        );
    }
    const int atomCount = molfileParserDetail::integerField(
        currLine, 0, 3, "atom count"
    );
    const int bondCount = molfileParserDetail::integerField(
        currLine, 3, 3, "bond count"
    );
    if (atomCount < 0 || bondCount < 0)
        throw std::runtime_error(
            "invalid molfile: atom and bond counts must be nonnegative"
        );

    if (verbose)
    {
        std::cout << "Molfile: " << atomCount << " atoms, " << bondCount
                  << " bonds\n";
    }
    parsed.atoms.reserve(static_cast<std::size_t>(atomCount));

    for (int atomIndex = 0; atomIndex < atomCount; atomIndex++)
    {
        molfileParserDetail::readLine(molfile, currLine, "atom line");
        const std::string_view atomField = molfileParserDetail::trimSpaces(
            molfileParserDetail::fixedField(
                currLine, 31, 3, "atom line"
            )
        );
        if (atomField.empty())
            throw std::runtime_error("invalid molfile: missing atom type");

        std::string atomType(atomField);
        parsed.addAtom(std::move(atomType));
    }
    for (int bondIndex = 0; bondIndex < bondCount; bondIndex++)
    {
        molfileParserDetail::readLine(molfile, currLine, "bond line");
        const int atomA = molfileParserDetail::integerField(
            currLine, 0, 3, "first bond atom"
        );
        const int atomB = molfileParserDetail::integerField(
            currLine, 3, 3, "second bond atom"
        );
        const int bondOrder = molfileParserDetail::integerField(
            currLine, 6, 3, "bond order"
        );
        if (atomA == atomB)
            throw std::runtime_error(
                "invalid molfile: self-loop bonds are not supported"
            );
        if (
            atomA < 1 || atomA > atomCount ||
            atomB < 1 || atomB > atomCount
        )
        {
            throw std::runtime_error(
                "invalid molfile: bond endpoint is outside the atom range"
            );
        }
        if (bondOrder < 0)
            throw std::runtime_error(
                "invalid molfile: bond order must be nonnegative"
            );
        if (bondOrder == 0)
            throw std::runtime_error(
                "invalid molfile: zero-order bonds are not supported"
            );
        parsed.addBond(atomA - 1, atomB - 1, bondOrder);
    }
    if (removeHydrogens)
    {
        for (size_t atomIndex = 0; atomIndex < parsed.atoms.size(); atomIndex++)
        {
            if (parsed.atomType(atomIndex) == "H") parsed.removeAtom(atomIndex);
        }
        parsed.removeAndCollapse();
    }
    if (verbose) parsed.printToCout();
    molecule = std::move(parsed);
}
