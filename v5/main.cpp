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
 */
void writeoutIntermediateMAs(const string &filename)
{
    ofstream ofs(filename);
    const int compensation = disjointCompensation ? disjointFragments : 1;
    for (size_t i = 0; i < intermediateMAs.size(); i++)
    {
        ofs << intermediateMAs[i].first << ' '
            << intermediateMAs[i].second - compensation + 1 << '\n';
    }
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
    const clock_t calculationStart = clock();
    outputFile << outputBase << " has assembly index: ";
    improvedBnB(mol_graph, outputFile);
    outputFile << "time to completion: " << clock() - calculationStart << '\n';
    if (writeIntermediateMAs) writeoutIntermediateMAs(outputBase + "IntermediateMAs");
    return true;
}

/**
 * @brief Memory usage tracker, works for linux only
 * 
 * @param outputFilename output filename
 */
void maxMemoryUsage(const string& outputFilename) {
    ifstream status_file("/proc/self/status");
    string line, peakMemory;

    while (getline(status_file, line))
    {
        if (line.rfind("VmPeak:", 0) == 0)
        {
            peakMemory = line;
            break;
        }
    }

    ofstream outFile(outputFilename);
    if (outFile.is_open()) {
        outFile << peakMemory << '\n';
        outFile.close();
    } 
    else
    {
        cerr << "Error: could not open output file.\n";
    }
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

    const bool succeeded = assemblyCalculator(arguments.input);

    #ifdef __linux__
        if (succeeded && memTest) maxMemoryUsage("memUsage");
    #endif

    return succeeded ? 0 : 1;
}
