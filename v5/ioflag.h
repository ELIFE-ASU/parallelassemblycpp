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
    removeHydrogens,
    compensateDisjoint,
    memoryReport,
    writeIntermediateMAs
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
            "Cooperatively stop the search after this many elapsed std::clock ticks.",
            {"runTime"}
        },
        {
            InputFlag::enumMax,
            "enum-max",
            "COUNT",
            "50000000",
            "Cap unique connected masks, including one-edge masks, in the initial DAG.",
            {"enumMax"}
        },
        {
            InputFlag::pathway,
            "pathway",
            "0|1",
            "1",
            "Write the recovered assembly pathway to INPUTPathway.",
            {}
        },
        {
            InputFlag::removeHydrogens,
            "remove-hydrogens",
            "0|1",
            "1",
            "Remove explicit hydrogen atoms from molfile inputs.",
            {"removeHydrogens"}
        },
        {
            InputFlag::compensateDisjoint,
            "compensate-disjoint",
            "0|1",
            "0",
            "Subtract one per processed-graph component after the first.",
            {"compensateDisjoint", "disjointCompensation"}
        },
        {
            InputFlag::memoryReport,
            "memory-report",
            "0|1",
            "0",
            "On Linux, write peak virtual memory (VmPeak) to ./memUsage.",
            {"memTest", "testMemory"}
        },
        {
            InputFlag::writeIntermediateMAs,
            "write-intermediate-mas",
            "0|1",
            "0",
            "Write each new best assembly index and its elapsed std::clock tick.",
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
            "option '" + inputFlagLabel(definition) +
            "' expects a non-negative whole number; received '" + value + "'"
        );
    }
    return parsed;
}

bool parseBooleanFlag(const InputFlagDefinition& definition, const string& value)
{
    if (value == "0") return false;
    if (value == "1") return true;
    throw std::invalid_argument(
        "option '" + inputFlagLabel(definition) +
        "' expects 0 or 1; received '" + value + "'"
    );
}

void applyInputFlag(const InputFlagDefinition& definition, const string& value)
{
    switch (definition.flag)
    {
        case InputFlag::runtime:
            runTimeMax = parseUnsignedFlag(definition, value);
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
                    "option '" + inputFlagLabel(definition) + "' expects a value from 1 to " +
                    std::to_string(std::numeric_limits<int>::max()) +
                    "; received '" + value + "'"
                );
            }
            ENUM_MAX = static_cast<int>(parsed);
            break;
        }

        case InputFlag::pathway:
            isPathway = parseBooleanFlag(definition, value);
            break;

        case InputFlag::removeHydrogens:
            removeHydrogens = parseBooleanFlag(definition, value);
            break;

        case InputFlag::compensateDisjoint:
            disjointCompensation = parseBooleanFlag(definition, value);
            break;

        case InputFlag::memoryReport:
            memTest = parseBooleanFlag(definition, value);
            break;

        case InputFlag::writeIntermediateMAs:
            writeIntermediateMAs = parseBooleanFlag(definition, value);
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

    for (int i = 1; i < argc; i++)
    {
        const string argument = argv[i];

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
                    "option '" + inputFlagLabel(*definition) + "' requires =<" +
                    definition->valueName + ">"
                );
            }

            const int flagId = static_cast<int>(definition->flag);
            if (!seenFlags.insert(flagId).second)
            {
                throw std::invalid_argument(
                    "option '" + inputFlagLabel(*definition) + "' was specified more than once"
                );
            }
            applyInputFlag(*definition, option.substr(equals + 1));
            continue;
        }

        if (!parsed.input.empty())
        {
            throw std::invalid_argument(
                "expected one INPUT, but received an extra argument '" + argument + "'"
            );
        }
        parsed.input = argument;
    }

    if (!parsed.showHelp && parsed.input.empty())
    {
        throw std::invalid_argument("missing required INPUT");
    }
    return parsed;
}
