#include "parallelassemblycpp.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}
}

int main()
{
    using namespace parallelassemblycpp;
    try
    {
        const std::vector<std::string> inputs{
            "", "x", "abababab", "abcddcba", "\xc3\xa9\xc3\xa9\xc3\xa9\xc3\xa9",
            std::string("a\0a\0", 4), "a\na\n"
        };
        for (bool reversed : {false, true})
        {
            StringCalculationOptions options;
            options.acceptReversed = reversed;
            const auto exact = calculateStringBatch(inputs, options);
            const auto exactBatch = calculateStringBatch({"abab", std::string("\xff", 1), "aaaa"}, options);
            require(exactBatch[0].succeeded && !exactBatch[1].succeeded && exactBatch[2].succeeded,
                    "exact batch did not recover after invalid UTF-8");
            options.rePairUpperBound = true;
            const auto bounds = calculateStringBatch(inputs, options);
            require(bounds.size() == inputs.size(), "batch size mismatch");
            for (std::size_t i = 0; i < bounds.size(); ++i)
            {
                require(exact[i].succeeded && bounds[i].succeeded, "string calculation failed");
                require(!exact[i].upperBoundOnly && bounds[i].upperBoundOnly, "wrong proof status");
                require(bounds[i].assemblyIndex >= exact[i].assemblyIndex, "bound below minimum");
                require(bounds[i].input == "<string>", "incorrect input marker");
                require(!bounds[i].runtimeLimitReached && !bounds[i].enumerationLimitReached,
                        "unexpected limit status");
            }
            require(bounds[0].assemblyIndex == -1, "incorrect empty bound");
            require(bounds[1].assemblyIndex == 0, "incorrect singleton bound");
            require(bounds[2].assemblyIndex == 3, "incorrect repeated block bound");
            require(bounds[4].assemblyIndex == 2, "UTF-8 bytes counted as symbols");
            const auto batch = calculateStringBatch({"abab", std::string("\xff", 1), "aaaa"}, options);
            require(batch[0].succeeded && !batch[1].succeeded && batch[2].succeeded,
                    "batch did not recover after invalid UTF-8");
            require(batch[1].error.find("UTF-8") != std::string::npos, "missing UTF-8 error");
            options.runtimeTicks = 0;
            const auto rejected = calculateString("abab", options);
            require(!rejected.succeeded && rejected.error.find("runtimeTicks") != std::string::npos,
                    "Re-Pair accepted a finite budget");
            options.rePairUpperBound = false;
            const auto limited = calculateString("abcabcabcabc", options);
            require(limited.succeeded && limited.runtimeLimitReached,
                    "exact string runtime limit was ignored");
            options.runtimeTicks = std::numeric_limits<std::uint64_t>::max();
            const auto recovered = calculateString("abababab", options);
            require(recovered.succeeded && !recovered.runtimeLimitReached &&
                    !recovered.upperBoundOnly && recovered.assemblyIndex == 3,
                    "limited call contaminated a later exact calculation");
        }
        require(calculateStringBatch({}).empty(), "empty batch failed");
        std::cout << "string library checks passed\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
