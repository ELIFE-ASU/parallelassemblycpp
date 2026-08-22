#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <csignal>
#include <cstdlib>
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
#if defined(ASSEMBLYCPP_USE_OPENMP)
    #include <omp.h>
#endif
#if defined(ASSEMBLYCPP_USE_MPI)
    #include <mpi.h>
#endif
#ifdef _WIN32
    #include <windows.h>
#endif

#if defined(ASSEMBLYCPP_USE_OPENMP)
    #define ASSEMBLYCPP_SEARCH_LOCAL thread_local
#else
    #define ASSEMBLYCPP_SEARCH_LOCAL
#endif

#if defined(ASSEMBLYCPP_USE_MPI)
int assemblyCppMpiRank = 0;
int assemblyCppMpiSize = 1;
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

#ifndef ASSEMBLYCPP_LIBRARY_BUILD

bool isPrimaryProcess()
{
#if defined(ASSEMBLYCPP_USE_MPI)
    return assemblyCppMpiRank == 0;
#else
    return true;
#endif
}

#if defined(ASSEMBLYCPP_USE_OPENMP) || defined(ASSEMBLYCPP_USE_MPI)

struct ParallelReplicaResult
{
    int assemblyIndex = std::numeric_limits<int>::max();
    bool started = false;
    bool succeeded = false;
    bool runtimeLimitReached = false;
    bool enumerationLimitReached = false;
    string error;
};

size_t parallelMinimumBonds()
{
    constexpr size_t defaultMinimum = 32;
    const char *configured = std::getenv("ASSEMBLYCPP_PARALLEL_MIN_BONDS");
    if (configured == nullptr || *configured == '\0') return defaultMinimum;

    size_t result = 0;
    const char *end = configured + std::char_traits<char>::length(configured);
    const auto parsed = std::from_chars(configured, end, result);
    if (parsed.ec != std::errc() || parsed.ptr != end) return defaultMinimum;
    return result;
}

int localParallelThreadCount()
{
#if defined(ASSEMBLYCPP_USE_OPENMP)
    omp_set_dynamic(0);
    return max(1, min(omp_get_max_threads(), omp_get_thread_limit()));
#else
    return 1;
#endif
}

bool hasMultipleParallelWorkers()
{
    const int localThreads = localParallelThreadCount();
#if defined(ASSEMBLYCPP_USE_MPI)
    return assemblyCppMpiSize > 1 || localThreads > 1;
#else
    return localThreads > 1;
#endif
}

bool parallelSearchEligible(const molGraph &graph)
{
    // Pathway and improvement-history output require one deterministic winning
    // search. Keep those modes serial until witness transfer is implemented.
    return
        !isPathway &&
        !writeIntermediateMAs &&
        runTimeMax == std::numeric_limits<unsigned long long>::max() &&
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        !searchTelemetryEnabled &&
#endif
        hasMultipleParallelWorkers() &&
        static_cast<size_t>(graph.totalBonds) >= parallelMinimumBonds();
}

#if defined(ASSEMBLYCPP_USE_MPI)
bool mpiCommandLineOptionsAgree(const CommandLineArguments &arguments)
{
    constexpr size_t optionCount = 10;
    const array<unsigned long long, optionCount> local = {
        arguments.showHelp ? 1ULL : 0ULL,
        static_cast<unsigned long long>(ENUM_MAX),
        runTimeMax,
        isPathway ? 1ULL : 0ULL,
        removeHydrogens ? 1ULL : 0ULL,
        verbose ? 1ULL : 0ULL,
        disjointCompensation ? 1ULL : 0ULL,
        memTest ? 1ULL : 0ULL,
        writeIntermediateMAs ? 1ULL : 0ULL,
#ifdef ASSEMBLY_ENABLE_TELEMETRY
        searchTelemetryEnabled ? 1ULL : 0ULL
#else
        0ULL
#endif
    };
    array<unsigned long long, optionCount> minimum = local;
    array<unsigned long long, optionCount> maximum = local;
    MPI_Allreduce(
        local.data(),
        minimum.data(),
        static_cast<int>(local.size()),
        MPI_UNSIGNED_LONG_LONG,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        local.data(),
        maximum.data(),
        static_cast<int>(local.size()),
        MPI_UNSIGNED_LONG_LONG,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    return minimum == maximum;
}

uint64_t parallelGraphFingerprint(const molGraph &graph)
{
    uint64_t result = UINT64_C(1469598103934665603);
    auto mix = [&](uint64_t value)
    {
        result ^= value;
        result *= UINT64_C(1099511628211);
    };
    mix(static_cast<uint64_t>(graph.totalBonds));
    mix(static_cast<uint64_t>(graph.mg.size()));
    for (const atom &entry : graph.mg)
    {
        mix(static_cast<uint64_t>(entry.type.size()));
        for (const unsigned char value : entry.type) mix(value);
        mix(static_cast<uint64_t>(entry.list.size()));
        for (const bond &edge : entry.list)
        {
            mix(static_cast<uint16_t>(edge.n));
            mix(static_cast<uint16_t>(edge.type));
        }
    }
    return result;
}

bool mpiGraphsAgree(const molGraph &graph)
{
    const uint64_t local = parallelGraphFingerprint(graph);
    uint64_t minimum = local;
    uint64_t maximum = local;
    MPI_Allreduce(
        &local,
        &minimum,
        1,
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &local,
        &maximum,
        1,
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    return minimum == maximum;
}
#endif

bool runParallelSearch(molGraph &graph, ofstream &output)
{
    const int localThreads = localParallelThreadCount();
    int shardOffset = 0;
    int totalShards = localThreads;
#if defined(ASSEMBLYCPP_USE_MPI)
    MPI_Exscan(
        &localThreads,
        &shardOffset,
        1,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    if (assemblyCppMpiRank == 0) shardOffset = 0;
    MPI_Allreduce(
        &localThreads,
        &totalShards,
        1,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD
    );
#endif

    std::atomic<int> processBest(std::numeric_limits<int>::max());
    vector<ParallelReplicaResult> replicas(
        static_cast<size_t>(localThreads)
    );

    auto runReplica = [&](int threadIndex)
    {
        ParallelReplicaResult &result = replicas[threadIndex];
        result.started = true;
        searchShardIndex = static_cast<size_t>(shardOffset + threadIndex);
        searchShardCount = static_cast<size_t>(totalShards);
        sharedAssemblyIndex = &processBest;
        suppressSearchOutput = true;
        try
        {
            molGraph workerGraph = graph;
            ofstream discardedOutput;
            result.succeeded = improvedBnB(workerGraph, discardedOutput);
            result.assemblyIndex = lastCalculatedAssemblyIndex;
            result.runtimeLimitReached = runtimeLimitReached;
            result.enumerationLimitReached = enumerationLimitReached;
        }
        catch (const std::exception &exception)
        {
            result.error = exception.what();
        }
        catch (...)
        {
            result.error = "unknown parallel worker failure";
        }
        sharedAssemblyIndex = nullptr;
        suppressSearchOutput = false;
        searchShardIndex = 0;
        searchShardCount = 1;
    };

#if defined(ASSEMBLYCPP_USE_OPENMP)
    #pragma omp parallel num_threads(localThreads)
    {
        runReplica(omp_get_thread_num());
    }
#else
    runReplica(0);
#endif

    int localAssemblyIndex = std::numeric_limits<int>::max();
    int localSucceeded = 1;
    int localRuntimeLimit = 0;
    int localEnumerationLimit = 0;
    for (const ParallelReplicaResult &replica : replicas)
    {
        localAssemblyIndex = min(localAssemblyIndex, replica.assemblyIndex);
        localSucceeded &= replica.started && replica.succeeded ? 1 : 0;
        localRuntimeLimit |= replica.runtimeLimitReached ? 1 : 0;
        localEnumerationLimit |= replica.enumerationLimitReached ? 1 : 0;
    }

    int globalAssemblyIndex = localAssemblyIndex;
    int globalSucceeded = localSucceeded;
    int globalRuntimeLimit = localRuntimeLimit;
    int globalEnumerationLimit = localEnumerationLimit;
    int globalUserInterrupt = receivedUserInterrupt() ? 1 : 0;
#if defined(ASSEMBLYCPP_USE_MPI)
    MPI_Allreduce(
        &localAssemblyIndex,
        &globalAssemblyIndex,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &localSucceeded,
        &globalSucceeded,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &localRuntimeLimit,
        &globalRuntimeLimit,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce(
        &localEnumerationLimit,
        &globalEnumerationLimit,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    const int localUserInterrupt = globalUserInterrupt;
    MPI_Allreduce(
        &localUserInterrupt,
        &globalUserInterrupt,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
#endif

    if (globalUserInterrupt != 0)
    {
#ifdef _WIN32
        userInterruptReceived.store(true);
#else
        userInterruptReceived = 1;
#endif
    }

    lastCalculatedAssemblyIndex = globalAssemblyIndex;
    runtimeLimitReached = globalRuntimeLimit != 0;
    enumerationLimitReached = globalEnumerationLimit != 0;
    if (isPrimaryProcess() && globalSucceeded != 0)
    {
        output << globalAssemblyIndex << '\n';
        if (runtimeLimitReached) output << "status: runtime limit reached\n";
        if (enumerationLimitReached)
            output << "status: enumeration limit reached\n";
    }
    if (isPrimaryProcess())
    {
        for (const ParallelReplicaResult &replica : replicas)
        {
            if (!replica.error.empty())
                cerr << "error: parallel worker failed: " << replica.error << '\n';
        }
        if (globalSucceeded == 0)
            cerr << "error: one or more parallel workers did not complete\n";
    }
    return globalSucceeded != 0;
}

#endif

bool runConfiguredSearch(molGraph &graph, ofstream &output)
{
#if defined(ASSEMBLYCPP_USE_MPI)
    if (!mpiGraphsAgree(graph))
    {
        if (isPrimaryProcess())
            cerr << "error: MPI ranks loaded different input graphs\n";
        return false;
    }
    int useParallelSearch = 0;
    if (isPrimaryProcess())
        useParallelSearch = parallelSearchEligible(graph) ? 1 : 0;
    MPI_Bcast(&useParallelSearch, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (useParallelSearch != 0) return runParallelSearch(graph, output);
#elif defined(ASSEMBLYCPP_USE_OPENMP)
    if (parallelSearchEligible(graph)) return runParallelSearch(graph, output);
#endif

#if defined(ASSEMBLYCPP_USE_MPI)
    int succeeded = 1;
    if (isPrimaryProcess())
    {
        try
        {
            succeeded = improvedBnB(graph, output) ? 1 : 0;
        }
        catch (const std::exception &exception)
        {
            cerr << "error: serial MPI-root search failed: "
                 << exception.what() << '\n';
            succeeded = 0;
        }
    }
    MPI_Bcast(&succeeded, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return succeeded != 0;
#else
    return improvedBnB(graph, output);
#endif
}

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
    searchCancellationFlag.store(false);
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    resetSearchTelemetry();
#endif
    const bool explicitMolfile = hasMolfileExtension(input);
    const string outputBase =
        explicitMolfile ? input.substr(0, input.size() - 4) : input;
    molGraph mol_graph;
    string inputError;
#if defined(ASSEMBLYCPP_USE_MPI)
    const bool configuredVerbose = verbose;
    if (!isPrimaryProcess()) verbose = false;
#endif
    int inputLoaded = loadMoleculeInput(input, mol_graph, inputError) ? 1 : 0;
#if defined(ASSEMBLYCPP_USE_MPI)
    verbose = configuredVerbose;
    int allInputsLoaded = inputLoaded;
    MPI_Allreduce(
        &inputLoaded,
        &allInputsLoaded,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    inputLoaded = allInputsLoaded;
#endif
    if (!inputLoaded)
    {
        if (isPrimaryProcess())
        {
            if (inputError.empty()) inputError = "input failed on another MPI rank";
            cerr << "error: " << inputError << '\n';
        }
        return false;
    }

    const string outputName = outputBase + "Out";
    ofstream outputFile;
    int outputReady = 1;
    if (isPrimaryProcess())
    {
        outputFile.open(outputName);
        if (!outputFile.is_open()) outputReady = 0;
    }
#if defined(ASSEMBLYCPP_USE_MPI)
    MPI_Bcast(&outputReady, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (!outputReady)
    {
        if (isPrimaryProcess())
            cerr << "error: could not open output file '" << outputName << "'\n";
        return false;
    }

    moleculeName = outputBase + "Pathway";
    if (isPrimaryProcess())
        outputFile << outputBase << " has assembly index: ";
    // improvedBnB propagates recoverPathway2's requested-output status.
    bool calculationSucceeded = false;
    try
    {
        calculationSucceeded = runConfiguredSearch(mol_graph, outputFile);
    }
    catch (const std::exception &exception)
    {
        if (isPrimaryProcess())
            cerr << "error: calculation failed for '" << input << "': "
                 << exception.what() << '\n';
        return false;
    }
    bool outputsSucceeded = calculationSucceeded;

    if (
        isPrimaryProcess() && writeIntermediateMAs &&
        !writeoutIntermediateMAs(outputBase + "IntermediateMAs")
    )
    {
        outputsSucceeded = false;
    }

    if (isPrimaryProcess() && receivedUserInterrupt())
    {
        cout << "status: interrupted by user\n";
        outputFile << "status: interrupted by user\n";
    }
    if (isPrimaryProcess())
        outputFile << "time elapsed: " << elapsedClockTicks() << '\n';

    if (isPrimaryProcess()) outputFile.close();
    if (isPrimaryProcess() && !outputFile)
    {
        cerr << "error: could not write output file '" << outputName << "'\n";
        outputsSucceeded = false;
    }
#ifdef ASSEMBLY_ENABLE_TELEMETRY
    if (
        isPrimaryProcess() && searchTelemetryEnabled &&
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
    searchCancellationFlag.store(false);
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
#if defined(ASSEMBLYCPP_USE_MPI)
class AssemblyCppMpiSession
{
    bool initialized = false;
    bool usable = false;

public:
    AssemblyCppMpiSession(int &argc, char **&argv)
    {
        int provided = MPI_THREAD_SINGLE;
#if defined(ASSEMBLYCPP_USE_OPENMP)
        constexpr int required = MPI_THREAD_FUNNELED;
#else
        constexpr int required = MPI_THREAD_SINGLE;
#endif
        if (MPI_Init_thread(&argc, &argv, required, &provided) != MPI_SUCCESS)
            return;
        initialized = true;
        MPI_Comm_rank(MPI_COMM_WORLD, &assemblyCppMpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &assemblyCppMpiSize);
        usable = provided >= required;
    }

    ~AssemblyCppMpiSession()
    {
        if (initialized) MPI_Finalize();
    }

    explicit operator bool() const noexcept {return usable;}
};
#endif

int main(int argc, char** argv)
{
    ios::sync_with_stdio(false);
#if defined(ASSEMBLYCPP_USE_MPI)
    AssemblyCppMpiSession mpiSession(argc, argv);
    if (!mpiSession)
    {
        if (assemblyCppMpiRank == 0)
            cerr << "error: MPI could not provide the required thread level\n";
        return 1;
    }
#endif
    #ifdef _WIN32
        static_cast<void>(SetConsoleCtrlHandler(CtrlHandler, TRUE));
    #else
        signal(SIGINT, signalHandler);
    #endif
    CommandLineArguments arguments;
    string argumentError;
    int argumentsValid = 1;
    try
    {
        arguments = parseCommandLine(argc, argv);
    }
    catch (const std::invalid_argument& error)
    {
        argumentError = error.what();
        argumentsValid = 0;
    }

#if defined(ASSEMBLYCPP_USE_MPI)
    int allArgumentsValid = argumentsValid;
    MPI_Allreduce(
        &argumentsValid,
        &allArgumentsValid,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    if (allArgumentsValid == 0)
    {
        if (isPrimaryProcess())
        {
            if (argumentError.empty())
                argumentError = "command line was invalid on another MPI rank";
            cerr << "error: " << argumentError << "\n"
                 << "Usage: AssemblyCpp INPUT [OPTIONS]\n"
                 << "Try 'AssemblyCpp --help' for more information.\n";
        }
        return 2;
    }
    if (!mpiCommandLineOptionsAgree(arguments))
    {
        if (isPrimaryProcess())
            cerr << "error: MPI ranks received different command-line options\n";
        return 2;
    }
#else
    if (argumentsValid == 0)
    {
        cerr << "error: " << argumentError << "\n"
             << "Usage: AssemblyCpp INPUT [OPTIONS]\n"
             << "Try 'AssemblyCpp --help' for more information.\n";
        return 2;
    }
#endif

    if (arguments.showHelp)
    {
        if (isPrimaryProcess()) help();
        return 0;
    }

    bool succeeded = assemblyCalculator(arguments.input);

    // All search and pathway output is complete; ignore new interrupts while
    // final process-level reports and stream buffers are finished.
    disableInterruptHandler();

    #ifdef __linux__
        if (isPrimaryProcess() && succeeded && memTest)
            succeeded = maxMemoryUsage("memUsage");
    #endif

    if (!succeeded) return 1;
    if (receivedUserInterrupt()) return 130;
    return 0;
}
#endif
