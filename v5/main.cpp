#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <sys/stat.h>
#include <iomanip>
#include <bitset>
#include <csignal>
#include <ctime>
#include <vector>
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif

using namespace std;
typedef vector<int> vi;
typedef vector<bool> vb;
typedef pair<int, int> pii;
#define BITSET_LENGTH 512
#define MAX_INT 2147483647
#define HASH_DEPTH_MAX 7
using standardBitset = bitset<BITSET_LENGTH>;

#include "globalPrimitives.h"
#include "ufds.h"
#include "molGraph.h"
#include "vf2.h"
#include "molfileParser.h"
#include "graphio.h"
#include "treeCanon.h"
#include "assemblyState.h"
#include "graphHashes.h"
#include "dagEnumeration.h"
#include "duplicateMatching.h"
#include "fragmentation.h"
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
    return filename.size() >= 4 && filename.compare(filename.size() - 4, 4, ".mol") == 0;
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
    const bool explicitMolfile = hasMolfileExtension(input);
    const string outputBase = explicitMolfile ? input.substr(0, input.size() - 4) : input;
    const string molfileName = explicitMolfile ? input : input + ".mol";
    molGraph mol_graph;
    ifstream molfile(molfileName);

    if (molfile.is_open())
    {
        cout << "opening " << molfileName << '\n';
        molfileParser(molfile, mol_graph);
    }
    else if (!explicitMolfile)
    {
        ifstream graphFile(input);
        if (graphFile.is_open())
        {
            cout << "opening " << input << '\n';
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
    const bool calculationSucceeded = improvedBnB(mol_graph, outputFile);
    bool outputsSucceeded = calculationSucceeded;

    if (
        writeIntermediateMAs &&
        !writeoutIntermediateMAs(outputBase + "IntermediateMAs")
    )
    {
        outputsSucceeded = false;
    }

    // Once all potentially long-running output generation is complete, stop
    // accepting new interrupts and record any interrupt already observed.
    disableInterruptHandler();
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
    return outputsSucceeded;
}

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

int main(int argc, char** argv)
{
    #ifdef _WIN32
        if (SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE)) {}
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

    #ifdef __linux__
        if (succeeded && memTest) succeeded = maxMemoryUsage("memUsage");
    #endif

    if (!succeeded) return 1;
    if (receivedUserInterrupt()) return 130;
    return 0;
}
