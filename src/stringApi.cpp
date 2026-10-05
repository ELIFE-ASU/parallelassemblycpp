#include "parallelassemblycpp.h"

#include <ctime>
#include <exception>
#include <type_traits>

#include "stringAssembly.h"
#include "stringRepair.h"

namespace parallelassemblycpp
{

CalculationResult calculateString(
    std::string_view input,
    const StringCalculationOptions &options
)
{
    CalculationResult result;
    result.input = "<string>";
    if (options.rePairUpperBound &&
        options.runtimeTicks != std::numeric_limits<std::uint64_t>::max())
    {
        result.error = "runtimeTicks must be unlimited for rePairUpperBound";
        return result;
    }
    try
    {
        const std::clock_t started = std::clock();
        if (options.rePairUpperBound)
        {
            const auto bound = detail::stringRepair::calculate(input, options.acceptReversed);
            result.assemblyIndex = bound.upperBound;
            result.upperBoundOnly = true;
        }
        else
        {
            detail::stringAssembly::Options exactOptions;
            exactOptions.runtimeTicks = options.runtimeTicks;
            exactOptions.acceptReversed = options.acceptReversed;
            exactOptions.reconstructPathway = false;
            const auto exact = detail::stringAssembly::calculate(std::string(input), exactOptions);
            result.assemblyIndex = exact.assemblyIndex;
            result.runtimeLimitReached = exact.runtimeLimitReached;
        }
        const std::clock_t finished = std::clock();
        if (started != static_cast<std::clock_t>(-1) && finished != static_cast<std::clock_t>(-1))
        {
            using ClockTicks = std::make_unsigned_t<std::clock_t>;
            result.clockTicks =
                static_cast<ClockTicks>(finished) - static_cast<ClockTicks>(started);
        }
        result.succeeded = true;
    }
    catch (const std::exception &exception)
    {
        result.error = exception.what();
    }
    return result;
}

std::vector<CalculationResult> calculateStringBatch(
    const std::vector<std::string> &inputs,
    const StringCalculationOptions &options
)
{
    std::vector<CalculationResult> results;
    results.reserve(inputs.size());
    for (const auto &input : inputs) results.push_back(calculateString(input, options));
    return results;
}

} // namespace parallelassemblycpp
