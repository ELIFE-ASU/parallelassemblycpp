#include "parallelassemblycpp.h"

#include <algorithm>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Batch exact calculations through the public API. The independent Python
// checker supplies expected indices and deliberately mixes graph domains and
// preprocessing options in one process to exercise per-calculation state.
int main()
{
    try
    {
        std::string line;
        while (std::getline(std::cin, line))
        {
            std::replace(line.begin(), line.end(), '|', '\n');
            line.push_back('\n');
            std::cout << '[';
            for (int variant = 0; variant < 4; ++variant)
            {
                parallelassemblycpp::CalculationOptions options;
                options.removeHydrogens = (variant & 2) != 0;
                options.compensateDisjoint = (variant & 1) != 0;
                std::istringstream input(line);
                const auto result = parallelassemblycpp::calculateGraph(
                    input, options
                );
                if (!result.succeeded || result.runtimeLimitReached ||
                    result.enumerationLimitReached || result.upperBoundOnly)
                {
                    throw std::runtime_error(
                        "exact calculation failed or was limited: " + result.error
                    );
                }
                if (variant != 0) std::cout << ',';
                std::cout << result.assemblyIndex;
            }
            std::cout << "]\n";
        }
        if (std::cin.bad()) throw std::runtime_error("failed to read graph records");
    }
    catch (const std::exception &error)
    {
        std::cerr << "graph addition-chain probe: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
