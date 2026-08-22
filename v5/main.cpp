#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <csignal>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <memory_resource>
#include <new>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#ifdef _WIN32
    #include <windows.h>
#endif

#ifdef ASSEMBLYCPP_LIBRARY_BUILD
#include "assemblycpp.h"
#endif

using namespace std;
using vi = vector<int>;
using vb = vector<bool>;
using pii = pair<int, int>;
#include "activeWordMask.h"

#ifdef ASSEMBLYCPP_LIBRARY_BUILD
namespace assemblycpp::detail
{
#endif

constexpr int ceilLog2(int value)
{
    return std::bit_width(static_cast<unsigned int>(value - 1));
}

#include "globalPrimitives.h"
#ifdef ASSEMBLY_ENABLE_TELEMETRY
#include "searchTelemetry.h"
#endif
#include "ufds.h"
#include "molGraph.h"
#include "molfileParser.h"
#include "graphio.h"
#include "treeCanon.h"
#include "cyclicCanon.h"
#include "assemblyState.h"
#include "graphHashes.h"
#include "dagEnumeration.h"
#include "duplicateMatching.h"
#include "fragmentation.h"
#include "assemblyTranspositionTable.h"
#include "pathwayGenerator.h"
#include "improvedBnB.h"
#include "signalHandler.h"
#include "ioflag.h"
#include "help.h"

/**
 * @brief Function to write out intermediate MAs before the calculation has terminated
 * 
 * @param filename output filename
 * @return true if the complete output was written successfully.
 */
bool writeoutIntermediateMAs(const string &filename)
{
    ofstream ofs(filename);
    if (!ofs.is_open())
    {
        cerr << "error: could not open output file '" << filename << "'\n";
        return false;
    }

    for (size_t i = 0; i < intermediateMAs.size(); i++)
    {
        ofs << intermediateMAs[i].first << ' '
            << compensateDisjointAssemblyIndex(intermediateMAs[i].second) << '\n';
    }

    ofs.close();
    if (!ofs)
    {
        cerr << "error: could not write output file '" << filename << "'\n";
        return false;
    }
    return true;
}

bool hasMolfileExtension(const string &filename)
{
    return filename.ends_with(".mol");
}

/** Load one supported input without creating any output files. */
#ifdef ASSEMBLYCPP_LIBRARY_BUILD
bool loadMoleculeInput(
    const string &input,
    molGraph &molGraphOutput,
    string &error
)
{
    const bool explicitMolfile = hasMolfileExtension(input);
    const string molfileName = explicitMolfile ? input : input + ".mol";
    ifstream molfile(molfileName);

    try
    {
        if (molfile.is_open())
        {
            if (verbose) cout << "opening " << molfileName << '\n';
            molfileParser(molfile, molGraphOutput);
            return true;
        }
        if (!explicitMolfile)
        {
            ifstream graphFile(input);
            if (graphFile.is_open())
            {
                if (verbose) cout << "opening " << input << '\n';
                graphio(graphFile, molGraphOutput);
                return true;
            }
            error = "input file not found: '" + input + "' (also tried '" +
                    molfileName + "')";
            return false;
        }
        error = "input file not found: '" + input + "'";
        return false;
    }
    catch (const std::exception &exception)
    {
        error = "could not parse '" + input + "': " + exception.what();
        return false;
    }
}
#endif

#ifndef ASSEMBLYCPP_LIBRARY_BUILD
/**
 * @brief Read a molfile or native graph and calculate its assembly index.
 *
 * Molfile paths may include or omit the .mol extension. Native graph paths are
 * used exactly as provided.
 *
 * @param input Input path supplied on the command line.
 * @return true if the input was read and the calculation output was written.
 */
bool assemblyCalculator(const string &input)
{
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    resetSearchTelemetry();
#endif
    const bool explicitMolfile = hasMolfileExtension(input);
    const string outputBase =
        explicitMolfile ? input.substr(0, input.size() - 4) : input;
    const string molfileName = explicitMolfile ? input : input + ".mol";
    molGraph mol_graph;
    ifstream molfile(molfileName);

    try
    {
        if (molfile.is_open())
        {
            if (verbose) cout << "opening " << molfileName << '\n';
            molfileParser(molfile, mol_graph);
        }
        else if (!explicitMolfile)
        {
            ifstream graphFile(input);
            if (graphFile.is_open())
            {
                if (verbose) cout << "opening " << input << '\n';
                graphio(graphFile, mol_graph);
            }
            else
            {
                cerr << "error: input file not found: '" << input
                     << "' (also tried '" << molfileName << "')\n";
                return false;
            }
        }
        else
        {
            cerr << "error: input file not found: '" << input << "'\n";
            return false;
        }
    }
    catch (const std::exception &exception)
    {
        cerr << "error: could not parse '" << input << "': "
             << exception.what() << '\n';
        return false;
    }

    const string outputName = outputBase + "Out";
    ofstream outputFile(outputName);
    if (!outputFile.is_open())
    {
        cerr << "error: could not open output file '" << outputName << "'\n";
        return false;
    }

    moleculeName = outputBase + "Pathway";
    outputFile << outputBase << " has assembly index: ";
    // improvedBnB propagates recoverPathway2's requested-output status.
    bool calculationSucceeded = false;
    try
    {
        calculationSucceeded = improvedBnB(mol_graph, outputFile);
    }
    catch (const std::exception &exception)
    {
        cerr << "error: calculation failed for '" << input << "': "
             << exception.what() << '\n';
        return false;
    }
    bool outputsSucceeded = calculationSucceeded;

    if (
        writeIntermediateMAs &&
        !writeoutIntermediateMAs(outputBase + "IntermediateMAs")
    )
    {
        outputsSucceeded = false;
    }

    if (receivedUserInterrupt())
    {
        cout << "status: interrupted by user\n";
        outputFile << "status: interrupted by user\n";
    }
    outputFile << "time elapsed: " << elapsedClockTicks() << '\n';

    outputFile.close();
    if (!outputFile)
    {
        cerr << "error: could not write output file '" << outputName << "'\n";
        outputsSucceeded = false;
    }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (
        searchTelemetryEnabled &&
        !writeSearchTelemetry(outputBase + "Telemetry.json")
    )
    {
        outputsSucceeded = false;
    }
#endif
    return outputsSucceeded;
}
#endif

#ifdef ASSEMBLYCPP_LIBRARY_BUILD

namespace
{

class LibraryOptionScope
{
    int previousEnumerationLimit = ENUM_MAX;
    unsigned long long previousRuntimeTicks = runTimeMax;
    bool previousPathway = isPathway;
    bool previousRemoveHydrogens = removeHydrogens;
    bool previousCompensateDisjoint = disjointCompensation;
    bool previousVerbose = verbose;
    bool previousMemoryReport = memTest;
    bool previousIntermediateMAs = writeIntermediateMAs;

public:
    explicit LibraryOptionScope(const assemblycpp::CalculationOptions &options)
    {
        ENUM_MAX = options.enumerationLimit;
        runTimeMax = options.runtimeTicks;
        isPathway = false;
        removeHydrogens = options.removeHydrogens;
        disjointCompensation = options.compensateDisjoint;
        verbose = options.verbose;
        memTest = false;
        writeIntermediateMAs = false;
    }

    ~LibraryOptionScope()
    {
        ENUM_MAX = previousEnumerationLimit;
        runTimeMax = previousRuntimeTicks;
        isPathway = previousPathway;
        removeHydrogens = previousRemoveHydrogens;
        disjointCompensation = previousCompensateDisjoint;
        verbose = previousVerbose;
        memTest = previousMemoryReport;
        writeIntermediateMAs = previousIntermediateMAs;
    }
};

void prepareInProcessCalculation()
{
    // A runtime budget uses the shared stop flag. Clear it before each item so
    // a limited calculation cannot stop the remainder of an in-process batch.
#ifdef _WIN32
    interruptFlag.store(false);
#else
    interruptFlag = 0;
#endif
}

assemblycpp::CalculationResult calculateLoadedMolecule(
    molGraph &graph,
    const string &input
)
{
    assemblycpp::CalculationResult result;
    result.input = input;
    prepareInProcessCalculation();
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    resetSearchTelemetry();
#endif
    // An unopened stream is a portable no-file sink for the legacy internal
    // search writer. The public API returns the same value directly below.
    ofstream discardedOutput;
    result.succeeded = improvedBnB(graph, discardedOutput);
    result.assemblyIndex = lastCalculatedAssemblyIndex;
    result.clockTicks = elapsedClockTicks();
    result.runtimeLimitReached = runtimeLimitReached;
    result.enumerationLimitReached = enumerationLimitReached;
    if (!result.succeeded)
        result.error = "the requested calculation output could not be produced";
    return result;
}

bool validLibraryOptions(
    const assemblycpp::CalculationOptions &options,
    string &error
)
{
    if (options.enumerationLimit >= 1) return true;
    error = "enumerationLimit must be at least one";
    return false;
}

} // namespace

} // namespace assemblycpp::detail
using namespace assemblycpp::detail;

assemblycpp::CalculationResult assemblycpp::calculateMolfile(
    std::istream &molfile,
    const CalculationOptions &options
)
{
    CalculationResult result;
    result.input = "<stream>";
    if (!validLibraryOptions(options, result.error)) return result;

    LibraryOptionScope optionScope(options);
    try
    {
        molGraph graph;
        molfileParser(molfile, graph);
        return calculateLoadedMolecule(graph, result.input);
    }
    catch (const std::exception &exception)
    {
        result.error = exception.what();
        return result;
    }
}

assemblycpp::CalculationResult assemblycpp::calculate(
    const std::string &input,
    const CalculationOptions &options
)
{
    CalculationResult result;
    result.input = input;
    if (!validLibraryOptions(options, result.error)) return result;

    LibraryOptionScope optionScope(options);
    try
    {
        molGraph graph;
        if (!loadMoleculeInput(input, graph, result.error)) return result;
        return calculateLoadedMolecule(graph, input);
    }
    catch (const std::exception &exception)
    {
        result.error = exception.what();
        return result;
    }
}

std::vector<assemblycpp::CalculationResult> assemblycpp::calculateBatch(
    const std::vector<std::string> &inputs,
    const CalculationOptions &options
)
{
    std::vector<CalculationResult> results;
    results.reserve(inputs.size());
    for (const string &input : inputs)
    {
        results.push_back(calculate(input, options));
    }
    return results;
}

#endif

#ifndef ASSEMBLYCPP_LIBRARY_BUILD
/**
 * @brief Memory usage tracker, works for linux only
 * 
 * @param outputFilename output filename
 * @return true if VmPeak was read and the complete report was written.
 */
bool maxMemoryUsage(const string& outputFilename)
{
    const string statusFilename = "/proc/self/status";
    ifstream status_file(statusFilename);
    if (!status_file.is_open())
    {
        cerr << "error: could not open memory status file '"
             << statusFilename << "'\n";
        return false;
    }

    string line, peakMemory;

    while (getline(status_file, line))
    {
        if (line.rfind("VmPeak:", 0) == 0)
        {
            peakMemory = line;
            break;
        }
    }

    if (status_file.bad())
    {
        cerr << "error: could not read memory status file '"
             << statusFilename << "'\n";
        return false;
    }
    if (peakMemory.empty())
    {
        cerr << "error: memory status file '" << statusFilename
             << "' does not contain 'VmPeak:'\n";
        return false;
    }

    ofstream outFile(outputFilename);
    if (!outFile.is_open())
    {
        cerr << "error: could not open output file '" << outputFilename << "'\n";
        return false;
    }

    outFile << peakMemory << '\n';
    outFile.close();
    if (!outFile)
    {
        cerr << "error: could not write output file '" << outputFilename << "'\n";
        return false;
    }
    return true;
}
#endif

#ifndef ASSEMBLYCPP_NO_MAIN
int main(int argc, char** argv)
{
    ios::sync_with_stdio(false);
    #ifdef _WIN32
        static_cast<void>(SetConsoleCtrlHandler(CtrlHandler, TRUE));
    #else
        signal(SIGINT, signalHandler);
    #endif
    CommandLineArguments arguments;
    try
    {
        arguments = parseCommandLine(argc, argv);
    }
    catch (const std::invalid_argument& error)
    {
        cerr << "error: " << error.what() << "\n"
             << "Usage: AssemblyCpp INPUT [OPTIONS]\n"
             << "Try 'AssemblyCpp --help' for more information.\n";
        return 2;
    }

    if (arguments.showHelp)
    {
        help();
        return 0;
    }

    bool succeeded = assemblyCalculator(arguments.input);

    // All search and pathway output is complete; ignore new interrupts while
    // final process-level reports and stream buffers are finished.
    disableInterruptHandler();

    #ifdef __linux__
        if (succeeded && memTest) succeeded = maxMemoryUsage("memUsage");
    #endif

    if (!succeeded) return 1;
    if (receivedUserInterrupt()) return 130;
    return 0;
}
#endif
