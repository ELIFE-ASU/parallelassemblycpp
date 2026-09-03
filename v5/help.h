#pragma once

/**
 * @brief Print command-line usage and option documentation.
 */
#ifdef ASSEMBLY_ENABLE_TELEMETRY
#define ASSEMBLY_TELEMETRY_OUTPUT_HELP \
    "  INPUTTelemetry.json   Search telemetry (--telemetry=1).\n"
#else
#define ASSEMBLY_TELEMETRY_OUTPUT_HELP ""
#endif

void help()
{
    cout << R"(AssemblyCpp v5

Usage:
  AssemblyCpp INPUT [OPTIONS]
  AssemblyCpp --help

Input:
  A V2000 MOL/SDF file or an AssemblyCpp native graph file.
  Existing .mol and .sdf suffixes are matched case-insensitively.
  A missing suffix tries the lowercase .mol spelling.
  An SDF input reads its first V2000 structure.
  Molfile output names omit a recognised suffix.

Options:
  -h, --help
      Show this help and exit.
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
  Boolean values are 0 or 1.
  In a parallel-enabled executable, --parallel=auto uses a work estimate from
  the prepared root jobs and DAG. If it selects serial execution, it reports
  the reason. --parallel=on bypasses the estimate, but fails when parallel
  execution cannot be honored. --parallel=off always runs serially.
  --threads sets the local thread count for each process; auto uses the OpenMP
  runtime default.
  Finite --runtime budgets and --write-intermediate-mas require serial search;
  a parallel-enabled executable reports an auto fallback or an on-mode error.
  Pathway output is supported after parallel optimization by deterministic
  reconstruction of a winning pathway.
  --runtime is a cooperative std::clock budget and may overrun while an
  operation finishes. CLOCKS_PER_SEC converts ticks to seconds; the clock
  source is platform-specific.
  --enum-max includes one-edge masks.
  A limited search records its best index and status in INPUTOut; the index may
  not be minimal.

Outputs:
  INPUTOut              Assembly index, status, and std::clock ticks.
  INPUTPathway          Recovered pathway JSON (--pathway=1).
  INPUTIntermediateMAs  Improved indices and ticks when enabled.
)" ASSEMBLY_TELEMETRY_OUTPUT_HELP R"(  ./memUsage            Linux VmPeak report (--memory-report=1).

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
  AssemblyCpp molecule.mol
  AssemblyCpp molecule --pathway=0 --enum-max=1000000
)";
}

#undef ASSEMBLY_TELEMETRY_OUTPUT_HELP
