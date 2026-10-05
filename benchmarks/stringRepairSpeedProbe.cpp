// One literal UTF-8 string per process; timing excludes file I/O and includes
// decoding. Use an external process timeout when comparing expensive searches.
#include "../src/stringAssembly.h"
#include "../src/stringRepair.h"

#include <chrono>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv)
{
    if ((argc != 3 && argc != 4) ||
        (std::string_view(argv[1]) != "exact" && std::string_view(argv[1]) != "re-pair") ||
        (argc == 4 && std::string_view(argv[3]) != "--accept-reversed"))
    {
        std::cerr << "usage: stringRepairSpeedProbe exact|re-pair INPUT [--accept-reversed]\n";
        return 2;
    }
    try
    {
        std::ifstream input(argv[2], std::ios::binary);
        if (!input) throw std::runtime_error("could not open input");
        const std::string value{std::istreambuf_iterator<char>(input), {}};
        if (input.bad()) throw std::runtime_error("could not read input");
        const bool rePair = std::string_view(argv[1]) == "re-pair";
        const bool reversed = argc == 4;
        const auto length = parallelassemblycpp::detail::stringAssembly::implementation::decodeInput(value).size();
        std::cout << std::setprecision(17) << std::boolalpha
                  << "{\"event\":\"started\",\"method\":\"" << argv[1]
                  << "\",\"length\":" << length
                  << ",\"accept_reversed\":" << reversed
                  << ",\"trivial_upper_bound\":" << static_cast<int>(length) - 1
                  << "}" << std::endl;
        const auto wallStart = std::chrono::steady_clock::now();
        const auto cpuStart = std::clock();
        int index = -1;
        int rules = 0;
        int residual = 0;
        if (rePair)
        {
            const auto bound = parallelassemblycpp::detail::stringRepair::calculate(value, reversed);
            index = bound.upperBound;
            rules = bound.ruleCount;
            residual = bound.remainingFragments;
        }
        else
        {
            parallelassemblycpp::detail::stringAssembly::Options options;
            options.acceptReversed = reversed;
            options.reconstructPathway = false;
            index = parallelassemblycpp::detail::stringAssembly::calculate(value, options).assemblyIndex;
        }
        const double cpu = static_cast<double>(std::clock() - cpuStart) / CLOCKS_PER_SEC;
        const double wall = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wallStart).count();
        std::cout << "{\"event\":\"finished\",\"method\":\"" << argv[1]
                  << "\",\"assembly_index\":" << index
                  << ",\"upper_bound_only\":" << rePair
                  << ",\"exact_completed\":" << !rePair
                  << ",\"algorithm_seconds\":" << wall
                  << ",\"cpu_seconds\":" << cpu;
        if (rePair)
            std::cout << ",\"rule_count\":" << rules
                      << ",\"remaining_fragments\":" << residual;
        std::cout << "}" << std::endl;
    }
    catch (const std::exception &error)
    {
        std::cerr << "string repair speed probe: " << error.what() << '\n';
        return 1;
    }
}
