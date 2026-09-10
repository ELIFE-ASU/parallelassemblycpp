#ifndef PARALLELASSEMBLYCPP_H
#define PARALLELASSEMBLYCPP_H

#include <cstdint>
#include <istream>
#include <limits>
#include <string>
#include <vector>

namespace parallelassemblycpp
{

#if \
    defined(PARALLELASSEMBLYCPP_LIBRARY_BUILD) && \
    defined(__GNUC__) && \
    !defined(__clang__)
#define PARALLELASSEMBLYCPP_PUBLIC __attribute__((externally_visible))
#else
#define PARALLELASSEMBLYCPP_PUBLIC
#endif

/** Options for one or more in-process assembly-index calculations. */
struct CalculationOptions
{
    std::uint64_t runtimeTicks = std::numeric_limits<std::uint64_t>::max();
    int enumerationLimit = 50000000;
    /** Remove explicit H vertices from MOL/SDF and native graph inputs. */
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
PARALLELASSEMBLYCPP_PUBLIC CalculationResult calculateMolfile(
    std::istream& molfile,
    const CalculationOptions& options = {}
);

/** Calculate directly from a native graph stream without creating files. */
PARALLELASSEMBLYCPP_PUBLIC CalculationResult calculateGraph(
    std::istream& graph,
    const CalculationOptions& options = {}
);

/** Read a native graph or first V2000 MOL/SDF structure without output files. */
PARALLELASSEMBLYCPP_PUBLIC CalculationResult calculate(
    const std::string& input,
    const CalculationOptions& options = {}
);

/**
 * Calculate several files sequentially in the current process.
 *
 * The implementation currently uses process-global search workspaces and is
 * therefore reusable but not thread-safe. One result is returned per input.
 */
PARALLELASSEMBLYCPP_PUBLIC std::vector<CalculationResult> calculateBatch(
    const std::vector<std::string>& inputs,
    const CalculationOptions& options = {}
);

#undef PARALLELASSEMBLYCPP_PUBLIC

} // namespace parallelassemblycpp

#endif
