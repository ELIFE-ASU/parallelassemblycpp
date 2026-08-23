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
  A V2000 molfile (.mol suffix optional) or an AssemblyCpp native graph file.
  Molfile output names omit the .mol suffix.

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
        for (size_t i = 0; i < definition.aliases.size(); i++)
        {
            if (i > 0) cout << ", ";
            cout << "-" << definition.aliases[i] << "=<" << definition.valueName << ">";
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
