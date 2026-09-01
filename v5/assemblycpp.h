#ifndef ASSEMBLYCPP_H
#define ASSEMBLYCPP_H

#include <cstdint>
#include <istream>
#include <limits>
#include <string>
#include <vector>

namespace assemblycpp
{

#if \
    defined(ASSEMBLYCPP_LIBRARY_BUILD) && \
    defined(__GNUC__) && \
    !defined(__clang__)
#define ASSEMBLYCPP_PUBLIC __attribute__((externally_visible))
#else
#define ASSEMBLYCPP_PUBLIC
#endif

/** Options for one or more in-process assembly-index calculations. */
struct CalculationOptions
{
    std::uint64_t runtimeTicks = std::numeric_limits<std::uint64_t>::max();
    int enumerationLimit = 50000000;
    bool removeHydrogens = true;
    bool compensateDisjoint = false;
    bool verbose = false;
};

/** Result returned without requiring callers to parse an output file. */
struct CalculationResult
{
    std::string input;
    int assemblyIndex = -1;
    std::uint64_t clockTicks = 0;
    bool succeeded = false;
    bool runtimeLimitReached = false;
    bool enumerationLimitReached = false;
    std::string error;

    explicit operator bool() const noexcept { return succeeded; }
};

/** Calculate directly from a V2000 molfile stream without creating files. */
ASSEMBLYCPP_PUBLIC CalculationResult calculateMolfile(
    std::istream& molfile,
    const CalculationOptions& options = {}
);

/** Read a molfile or native graph and calculate without creating output files. */
ASSEMBLYCPP_PUBLIC CalculationResult calculate(
    const std::string& input,
    const CalculationOptions& options = {}
);

/**
 * Calculate several files sequentially in the current process.
 *
 * The implementation currently uses process-global search workspaces and is
 * therefore reusable but not thread-safe. One result is returned per input.
 */
ASSEMBLYCPP_PUBLIC std::vector<CalculationResult> calculateBatch(
    const std::vector<std::string>& inputs,
    const CalculationOptions& options = {}
);

#undef ASSEMBLYCPP_PUBLIC

} // namespace assemblycpp

#endif
