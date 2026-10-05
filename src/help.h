#pragma once

#ifdef ASSEMBLY_ENABLE_TELEMETRY
#define ASSEMBLY_TELEMETRY_OUTPUT_HELP \
    "  INPUTTelemetry.json   Search telemetry (--telemetry=1).\n"
#else
#define ASSEMBLY_TELEMETRY_OUTPUT_HELP ""
#endif

/** Print command-line usage and option documentation for this build. */
void help()
{
    cout << R"(ParallelAssemblyCpp

Usage:
  ParallelAssemblyCpp INPUT [OPTIONS]
  ParallelAssemblyCpp [OPTIONS] -- INPUT
  ParallelAssemblyCpp --help

Input:
  By default, a V2000 MOL/SDF file or a ParallelAssemblyCpp native graph file.
  Case-insensitive .mol and .sdf suffixes select V2000 parsing. Otherwise an
  existing path is read as a native graph; a missing path is retried with .mol
  appended.
  An SDF input reads only its first record, which must be V2000.
  V2000 graph labels use atom symbols and bond orders; coordinates, charges,
  isotopes, stereochemistry, and other property fields are ignored.
  Molfile output names omit a recognised suffix.
  With --run-strings=1, INPUT is read exactly as a text file containing one
  UTF-8 string per line. Each Unicode code point is one symbol, without
  normalization. LF, CRLF, empty lines, and a final line without a newline
  are supported. String pathway positions and lengths count code points.

Options:
  -h, --help
      Show this help and exit.
  --
      Stop parsing options. The next argument is INPUT; extra inputs are errors.
)";

    for (const InputFlagDefinition& definition : inputFlagDefinitions())
    {
        cout << "  --" << definition.name << "=<" << definition.valueName << ">\n"
             << "      " << definition.description
             << " Default: " << definition.defaultValue << ".\n";
    }

    cout << R"(
Notes:
  Options may appear before or after INPUT. Use --name=value.
  Boolean values are 0 or 1. Each option may be specified only once, including
  when using an alias.
  --algorithm=full runs exact search (the default). --algorithm=re-pair returns
  a graph or string upper bound and does not prove the minimum. Re-Pair uses
  serial execution: --parallel=auto reports a fallback and --parallel=on fails.
  Re-Pair rejects explicit --runtime or --enum-max and enabled telemetry or
  intermediate-index output. Its pathway JSON is a replayable bound certificate.
  --upper-bound=graph-repair remains graph-only; do not combine it with --algorithm.
  A lone -- ends option parsing, so an INPUT whose name begins with a dash
  must follow it. Nothing after -- is an option: a later --help or
  --name=value is read as INPUT instead.
  In a parallel-enabled executable, --parallel=auto uses a work estimate from
  the prepared root jobs and DAG for molecular search. If it selects serial
  execution, it reports the reason. --parallel=on forces parallel search, but
  fails when parallel execution cannot be honored. --parallel=off always runs
  serially. --threads sets the local thread count for each process. In molecular
  search, auto caps the OpenMP runtime default by estimated work per worker,
  shared across MPI ranks with at least one thread per rank (at least two
  total workers for --parallel=on).
  Explicit counts are never reduced by the workload cap.
  --threads is unused when --parallel=off.
  Finite --runtime budgets and --write-intermediate-mas require serial search;
  a parallel-enabled executable reports an auto fallback or an on-mode error.
  Pathway output is supported after parallel optimization by deterministic
  reconstruction of a winning pathway.
  --runtime is a cooperative std::clock budget and may overrun while an
  operation finishes. CLOCKS_PER_SEC converts ticks to seconds; the clock
  source is platform-specific. In string mode the budget applies to each line.
  --enum-max includes one-edge masks and applies only to graph inputs.
  A limited search records its best index and status in INPUTOut; the index may
  not be minimal.
  Full string assembly distributes search branches within each line across
  OpenMP threads and/or MPI ranks, preserving input order and deterministic pathways.
  String threads are capped by available root jobs; short strings may not
  benefit from parallel execution. Telemetry and intermediate-index output
  are unavailable in string mode. Explicit --enum-max and --remove-hydrogens
  options, and --compensate-disjoint=1, are also rejected in string mode.
  --accept-palindromes=1 identifies a fragment with its reversal in string
  mode and is rejected for graph inputs. Disabled boolean options such as
  --compensate-disjoint=0 are accepted in either mode.
  --memory-report=1 is available only on Linux and rejects an input that
  would be overwritten by ./memUsage.

Outputs:
  INPUTOut              Index or upper bound, status when needed, and clock ticks.
  INPUTPathway          Graph pathway JSON (--pathway=1).
  INPUT_N_Pathway       String pathway JSON for zero-based line N (--pathway=1).
  INPUTIntermediateMAs  Elapsed ticks then improved index, one pair per line.
)" ASSEMBLY_TELEMETRY_OUTPUT_HELP R"(  ./memUsage            Linux VmPeak report (--memory-report=1).

Exit status:
  0    Help shown or requested outputs written, including limited searches.
  1    Input, calculation, or output failure.
  2    Invalid command line.
  130  Interrupted by the user after available outputs were written.
  Inspect INPUTOut for limits or a heuristic status before assuming minimality.

Legacy options:
  Canonical and legacy names accept one or two leading dashes.
  Renamed legacy options:
)";

    for (const InputFlagDefinition& definition : inputFlagDefinitions())
    {
        if (definition.aliases.empty()) continue;
        cout << "  --" << definition.name << ": ";
        for (size_t aliasIndex = 0; aliasIndex < definition.aliases.size(); aliasIndex++)
        {
            if (aliasIndex > 0) cout << ", ";
            cout << "-" << definition.aliases[aliasIndex] << "=<"
                << definition.valueName << ">";
        }
        cout << '\n';
    }

    cout << R"(
Examples:
  ParallelAssemblyCpp molecule.mol
  ParallelAssemblyCpp molecule.mol --algorithm=full
  ParallelAssemblyCpp molecule.mol --algorithm=re-pair
  ParallelAssemblyCpp molecule --pathway=0 --enum-max=1000000
  ParallelAssemblyCpp strings.txt --run-strings=1
  ParallelAssemblyCpp strings.txt --run-strings=1 --algorithm=re-pair
  ParallelAssemblyCpp --pathway=0 -- -dashed-name.mol
)";
}

#undef ASSEMBLY_TELEMETRY_OUTPUT_HELP
