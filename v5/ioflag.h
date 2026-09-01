#include <charconv>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

enum class InputFlag
{
    runtime,
    enumMax,
    pathway,
    parallel,
    threads,
    removeHydrogensFlag,
    verboseFlag,
    compensateDisjoint,
    memoryReport,
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    telemetry,
#endif
    writeIntermediateMAsFlag
};

struct InputFlagDefinition
{
    InputFlag flag;
    string name;
    string valueName;
    string defaultValue;
    string description;
    vector<string> aliases;
};

struct CommandLineArguments
{
    string input;
    bool showHelp = false;
};

/**
 * @brief Definitions shared by the parser and help output.
 */
const vector<InputFlagDefinition>& inputFlagDefinitions()
{
    static const vector<InputFlagDefinition> definitions = {
        {
            InputFlag::runtime,
            "runtime",
            "TICKS",
            "unlimited",
            "Stop after this many std::clock ticks.",
            {"runTime"}
        },
        {
            InputFlag::enumMax,
            "enum-max",
            "COUNT",
            "50000000",
            "Limit connected subgraphs in initial enumeration.",
            {"enumMax"}
        },
        {
            InputFlag::pathway,
            "pathway",
            "0|1",
            "1",
            "Write recovered pathway JSON to INPUTPathway.",
            {}
        },
        {
            InputFlag::parallel,
            "parallel",
            "auto|on|off",
            "off",
            "Select parallel search automatically, require it, or disable it.",
            {}
        },
        {
            InputFlag::threads,
            "threads",
            "auto|N",
            "auto",
            "Set the OpenMP thread count per process (N must be positive).",
            {}
        },
        {
            InputFlag::removeHydrogensFlag,
            "remove-hydrogens",
            "0|1",
            "1",
            "Remove explicit hydrogens from molfile inputs.",
            {"removeHydrogens"}
        },
        {
            InputFlag::verboseFlag,
            "verbose",
            "0|1",
            "0",
            "Print the parsed input graph.",
            {}
        },
        {
            InputFlag::compensateDisjoint,
            "compensate-disjoint",
            "0|1",
            "0",
            "Subtract one for each processed component after the first.",
            {"compensateDisjoint", "disjointCompensation"}
        },
        {
            InputFlag::memoryReport,
            "memory-report",
            "0|1",
            "0",
            "On Linux, write VmPeak to ./memUsage.",
            {"memTest", "testMemory"}
        },
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        {
            InputFlag::telemetry,
            "telemetry",
            "0|1",
            "0",
            "Write search telemetry to INPUTTelemetry.json.",
            {}
        },
#endif
        {
            InputFlag::writeIntermediateMAsFlag,
            "write-intermediate-mas",
            "0|1",
            "0",
            "Write index improvements to INPUTIntermediateMAs.",
            {"writeIntermediateMAs"}
        }
    };
    return definitions;
}

const InputFlagDefinition* findInputFlagDefinition(const string& name)
{
    for (const InputFlagDefinition& definition : inputFlagDefinitions())
    {
        if (name == definition.name) return &definition;
        for (const string& alias : definition.aliases)
        {
            if (name == alias) return &definition;
        }
    }
    return nullptr;
}

string inputFlagLabel(const InputFlagDefinition& definition)
{
    return "--" + definition.name;
}

unsigned long long parseUnsignedFlag(
    const InputFlagDefinition& definition,
    const string& value
)
{
    unsigned long long parsed = 0;
    const char* begin = value.data();
    const char* end = begin + value.size();
    const auto result = std::from_chars(begin, end, parsed);

    if (value.empty() || result.ec != std::errc() || result.ptr != end)
    {
        throw std::invalid_argument(
            inputFlagLabel(definition) +
            ": expected a non-negative integer; got '" + value + "'"
        );
    }
    return parsed;
}

bool parseBooleanFlag(const InputFlagDefinition& definition, const string& value)
{
    if (value == "0") return false;
    if (value == "1") return true;
    throw std::invalid_argument(
        inputFlagLabel(definition) + ": expected 0 or 1; got '" + value + "'"
    );
}

parallelMode parseParallelModeFlag(
    const InputFlagDefinition& definition,
    const string& value
)
{
    if (value == "auto") return parallelMode::automatic;
    if (value == "on") return parallelMode::on;
    if (value == "off") return parallelMode::off;
    throw std::invalid_argument(
        inputFlagLabel(definition) +
        ": expected auto, on, or off; got '" + value + "'"
    );
}

size_t parseParallelThreadCountFlag(
    const InputFlagDefinition& definition,
    const string& value
)
{
    if (value == "auto") return 0;

    const unsigned long long parsed = parseUnsignedFlag(definition, value);
    if (
        parsed < 1 ||
        parsed > static_cast<unsigned long long>(std::numeric_limits<int>::max())
    )
    {
        throw std::invalid_argument(
            inputFlagLabel(definition) +
            ": expected auto or an integer from 1 to " +
            std::to_string(std::numeric_limits<int>::max()) +
            "; got '" + value + "'"
        );
    }
    return static_cast<size_t>(parsed);
}

void applyInputFlag(const InputFlagDefinition& definition, const string& value)
{
    switch (definition.flag)
    {
        case InputFlag::runtime:
            maximumRuntimeTicks = parseUnsignedFlag(definition, value);
            break;

        case InputFlag::enumMax:
        {
            const unsigned long long parsed = parseUnsignedFlag(definition, value);
            if (
                parsed < 1 ||
                parsed > static_cast<unsigned long long>(std::numeric_limits<int>::max())
            )
            {
                throw std::invalid_argument(
                    inputFlagLabel(definition) + ": expected an integer from 1 to " +
                    std::to_string(std::numeric_limits<int>::max()) +
                    "; got '" + value + "'"
                );
            }
            maximumEnumerationCount = static_cast<int>(parsed);
            break;
        }

        case InputFlag::pathway:
            pathwayOutputEnabled = parseBooleanFlag(definition, value);
            break;

        case InputFlag::parallel:
            parallelExecutionMode = parseParallelModeFlag(definition, value);
            break;

        case InputFlag::threads:
            parallelThreadCount = parseParallelThreadCountFlag(definition, value);
            break;

        case InputFlag::removeHydrogensFlag:
            removeHydrogens = parseBooleanFlag(definition, value);
            break;

        case InputFlag::verboseFlag:
            verbose = parseBooleanFlag(definition, value);
            break;

        case InputFlag::compensateDisjoint:
            disjointCompensation = parseBooleanFlag(definition, value);
            break;

        case InputFlag::memoryReport:
            memoryReportEnabled = parseBooleanFlag(definition, value);
            break;

#ifdef ASSEMBLY_ENABLE_TELEMETRY
        case InputFlag::telemetry:
        {
            searchTelemetryEnabled = parseBooleanFlag(definition, value);
            break;
        }
#endif

        case InputFlag::writeIntermediateMAsFlag:
            writeIntermediateAssemblyIndices = parseBooleanFlag(definition, value);
            break;
    }
}

/**
 * @brief Parse one input path and any supported options.
 *
 * Options may appear before or after the input. Canonical options use
 * --kebab-case=value; former camelCase spellings and either one or two leading
 * dashes remain accepted for compatibility.
 *
 * @throws invalid_argument if the command line is incomplete or invalid.
 */
CommandLineArguments parseCommandLine(int argc, char** argv)
{
    CommandLineArguments parsed;
    std::unordered_set<int> seenFlags;
    bool optionsEnabled = true;

    for (int argumentIndex = 1; argumentIndex < argc; argumentIndex++)
    {
        const string argument = argv[argumentIndex];

        if (optionsEnabled && argument == "--")
        {
            optionsEnabled = false;
            continue;
        }
        if (optionsEnabled && (argument == "--help" || argument == "-h"))
        {
            parsed.showHelp = true;
            continue;
        }

        const bool isOption = optionsEnabled && argument.size() > 1 && argument[0] == '-';
        if (isOption)
        {
            const std::size_t prefixLength = argument.rfind("--", 0) == 0 ? 2 : 1;
            const string option = argument.substr(prefixLength);
            const std::size_t equals = option.find('=');
            const string name = option.substr(0, equals);
            const InputFlagDefinition* definition = findInputFlagDefinition(name);

            if (definition == nullptr)
            {
                throw std::invalid_argument("unknown option '" + argument + "'");
            }
            if (equals == string::npos)
            {
                throw std::invalid_argument(
                    inputFlagLabel(*definition) +
                    " requires a value in the form " +
                    inputFlagLabel(*definition) + "=<" +
                    definition->valueName + ">"
                );
            }

            const int flagId = static_cast<int>(definition->flag);
            if (!seenFlags.insert(flagId).second)
            {
                throw std::invalid_argument(
                    inputFlagLabel(*definition) + " may be specified only once"
                );
            }
            applyInputFlag(*definition, option.substr(equals + 1));
            continue;
        }

        if (!parsed.input.empty())
        {
            throw std::invalid_argument(
                "expected one INPUT; got extra argument '" + argument + "'"
            );
        }
        parsed.input = argument;
    }

    if (!parsed.showHelp && parsed.input.empty())
    {
        throw std::invalid_argument("INPUT is required");
    }
    return parsed;
}
