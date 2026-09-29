#include <iostream>
#include <sstream>
#include <string>
#include "parallelassemblycpp.h"

// Each stdin record is a five-line native graph with '|' replacing newlines.
int main()
{
    std::string line;
    while (std::getline(std::cin, line))
    {
        for (char &c : line) if (c == '|') c = '\n';
        std::istringstream input(line);
        parallelassemblycpp::CalculationOptions options;
        options.removeHydrogens = false;
        auto result = parallelassemblycpp::calculateGraph(input, options);
        for (char &c : result.error)
            if (c == '\n' || c == '\r' || c == '\t') c = ' ';
        std::cout << result.assemblyIndex << '\t' << result.succeeded << '\t'
                  << result.runtimeLimitReached << '\t'
                  << result.enumerationLimitReached << '\t' << result.error
                  << std::endl;
    }
    return std::cin.bad() ? 1 : 0;
}
