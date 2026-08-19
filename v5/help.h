/**
 * @brief Print command-line usage and option documentation.
 */
#ifdef ASSEMBLY_ENABLE_TELEMETRY
#define ASSEMBLY_TELEMETRY_OUTPUT_HELP \
    "  INPUTTelemetry.json   Search counters, cache rates, and phase memory when\n" \
    "                        --telemetry=1.\n"
#else
#define ASSEMBLY_TELEMETRY_OUTPUT_HELP ""
#endif

void help()
{
    cout << R"(AssemblyCpp v5

Usage:
  AssemblyCpp INPUT [OPTIONS]
  AssemblyCpp [OPTIONS] INPUT
  AssemblyCpp --help

Input:
  INPUT may be a .mol file, a molfile path with the .mol suffix omitted, or a
  file in AssemblyCpp's native graph format. For a molfile, output names use
  the input path without the .mol suffix.

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
Option notes:
  Boolean values must be exactly 0 (disabled) or 1 (enabled).
  --runtime is a cooperatively checked search budget measured in raw std::clock
  ticks. CLOCKS_PER_SEC ticks represent one second according to the platform's
  C++ runtime; whether clock() measures processor or elapsed time is
  implementation-specific. Checks occur between search operations, so the
  requested budget can be exceeded by the duration of an in-progress operation.
  A runtime or enumeration limit can stop an exhaustive search. INPUTOut still
  begins with the numeric best index found so far and then records which limit
  was reached; that value is not necessarily a proven minimum.

Outputs:
  INPUTOut              Assembly index, search status, and std::clock ticks.
  INPUTPathway          Recovered pathway when --pathway=1.
  INPUTIntermediateMAs  Improved intermediate indices when enabled.
)" ASSEMBLY_TELEMETRY_OUTPUT_HELP R"(  ./memUsage            Linux VmPeak report after other outputs succeed when
                        --memory-report=1.

Compatibility:
  Both one and two leading dashes are accepted for the canonical names and the
  former option names below:
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
  AssemblyCpp --remove-hydrogens=0 molecule.mol
)";
}

#undef ASSEMBLY_TELEMETRY_OUTPUT_HELP
