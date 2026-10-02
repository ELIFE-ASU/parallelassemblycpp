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

/**
 * Thread safety: the calculation entry points below read and write process-global
 * option and search state, so the API is reusable but not thread-safe. Only
 * one calculation may run at a time in a process; use separate processes for
 * concurrent work.
 */

/** Options for one or more in-process assembly-index calculations. */
struct CalculationOptions
{
    /** Cooperative std::clock budget for search; the maximum value is unlimited. */
    std::uint64_t runtimeTicks = std::numeric_limits<std::uint64_t>::max();
    /** Initial DAG mask limit, including one-edge masks; must be at least one. */
    int enumerationLimit = 50000000;
    /** Remove explicit H vertices from MOL/SDF and native graph inputs. */
    bool removeHydrogens = true;
    /** Subtract joins between processed disconnected components. */
    bool compensateDisjoint = false;
    /** Print parsing and search diagnostics to standard output. */
    bool verbose = false;
    /**
     * Return a GraphRePair-inspired heuristic upper bound instead of exact search.
     * Requires unlimited runtimeTicks; enumerationLimit is validated but unused.
     */
    bool graphRepairUpperBound = false;
};

/**
 * Result returned without requiring callers to parse an output file.
 *
 * succeeded reports whether a result was obtained, not whether minimality was
 * proved. A successful result is a proven minimum only when both limit flags
 * and upperBoundOnly are false.
 */
struct CalculationResult
{
    /** Supplied file name, or "<stream>" for either stream entry point. */
    std::string input;
    /** Best index found when succeeded is true. */
    int assemblyIndex = -1;
    /** Search duration in std::clock ticks, excluding input parsing. */
    std::uint64_t clockTicks = 0;
    bool succeeded = false;
    bool runtimeLimitReached = false;
    bool enumerationLimitReached = false;
    /** The result is a heuristic upper bound; minimality was not checked. */
    bool upperBoundOnly = false;
    /** Diagnostic for invalid options, input failures, or calculation failures. */
    std::string error;

    explicit operator bool() const noexcept { return succeeded; }
};

/**
 * Calculate from one V2000 header and atom/bond blocks without creating files.
 * Trailing molfile properties and subsequent SDF records are left unread.
 * Parsing, option, and calculation failures are reported in the result.
 */
PARALLELASSEMBLYCPP_PUBLIC CalculationResult calculateMolfile(
    std::istream& molfile,
    const CalculationOptions& options = {}
);

/**
 * Calculate from five native graph lines without creating files.
 * Trailing lines are left unread; failures are reported in the result.
 */
PARALLELASSEMBLYCPP_PUBLIC CalculationResult calculateGraph(
    std::istream& graph,
    const CalculationOptions& options = {}
);

/**
 * Read a native graph or MOL/SDF file without creating output files.
 *
 * Case-insensitive .mol and .sdf suffixes select V2000 parsing. Otherwise an
 * existing path is read as a native graph; a missing path is retried with .mol
 * appended. An SDF input reads only its first record, which must be V2000.
 * Parsing, option, and calculation failures are reported in the result.
 */
PARALLELASSEMBLYCPP_PUBLIC CalculationResult calculate(
    const std::string& input,
    const CalculationOptions& options = {}
);

/**
 * Calculate several files sequentially in the current process.
 *
 * One result is returned per input in the original order. A failed or limited
 * calculation does not prevent subsequent inputs from being processed.
 */
PARALLELASSEMBLYCPP_PUBLIC std::vector<CalculationResult> calculateBatch(
    const std::vector<std::string>& inputs,
    const CalculationOptions& options = {}
);

#undef PARALLELASSEMBLYCPP_PUBLIC

} // namespace parallelassemblycpp

#endif
