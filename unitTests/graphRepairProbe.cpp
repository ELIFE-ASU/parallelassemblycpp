// Emit graph-RePair certificates for the independent Python checker.
#define PARALLELASSEMBLYCPP_NO_MAIN
#include "../src/main.cpp"
#include "../src/graphRepair.h"

int main(int argc, char **argv)
{
    bool files = false;
    bool compensate = false;
    for (int index = 1; index < argc; ++index)
    {
        const string argument = argv[index];
        if (argument == "--files") files = true;
        else if (argument == "--compensate-disjoint") compensate = true;
        else
        {
            cerr << "unknown probe argument: " << argument << '\n';
            return 2;
        }
    }
    removeHydrogens = files;
    verbose = false;
    string line;
    while (getline(cin, line))
    {
        try
        {
            molGraph molecule;
            if (files)
            {
                string error;
                if (!loadMoleculeInput(line, molecule, error))
                    throw runtime_error(error);
            }
            else
            {
                replace(line.begin(), line.end(), '|', '\n');
                line.push_back('\n');
                istringstream input(line);
                graphio(input, molecule);
            }
            const auto started = chrono::steady_clock::now();
            const auto result = graphRepair::calculate(molecule, compensate);
            const double seconds = chrono::duration<double>(
                chrono::steady_clock::now() - started
            ).count();
            ostringstream certificate;
            graphRepair::writeJson(result, molecule, certificate);
            string document = certificate.str();
            const size_t closingBrace = document.find_last_not_of(" \t\n\r");
            if (closingBrace == string::npos || document[closingBrace] != '}')
                throw runtime_error("certificate is not a JSON object");
            cout << document.substr(0, closingBrace) << ",\"input\":";
            printJsonString(files ? line : "native graph", cout);
            cout << ",\"elapsed_seconds\":" << seconds << '}' << endl;
        }
        catch (const exception &error)
        {
            cerr << "graph repair probe: " << error.what() << '\n';
            return 1;
        }
    }
    return cin.bad() ? 1 : 0;
}
