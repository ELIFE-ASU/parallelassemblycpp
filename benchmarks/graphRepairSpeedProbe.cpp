// One calculation per process: graph_repair_speed.py times the whole process.
// The timers below exclude parsing and cover preparation plus calculation.
#define PARALLELASSEMBLYCPP_NO_MAIN
#include "../src/main.cpp"
#include <iomanip>

int main(int argc, char **argv)
{
    if (argc != 3 || (string(argv[1]) != "exact" && string(argv[1]) != "graph"))
    {
        cerr << "usage: graphRepairSpeedProbe exact|graph INPUT\n";
        return 2;
    }
    const string method = argv[1];
    verbose = false;
    suppressSearchOutput = true;
    pathwayOutputEnabled = false;
    removeHydrogens = true;
    disjointCompensation = false;
    maximumRuntimeTicks = numeric_limits<unsigned long long>::max();

    try
    {
        molGraph molecule;
        string error;
        if (!loadMoleculeInput(argv[2], molecule, error))
            throw runtime_error(error);
        // A flushed start event establishes that parsing finished, including
        // when the external process deadline terminates an unfinished search.
        cout << setprecision(17) << boolalpha
             << "{\"event\":\"started\",\"method\":\"" << method
             << "\",\"atoms\":" << molecule.atoms.size()
             << ",\"bonds\":" << molecule.totalBonds
             << ",\"components\":" << molecule.disjointFragments()
             << ",\"trivial_upper_bound\":" << max(0, molecule.totalBonds - 1)
             << ",\"enumeration_limit\":" << maximumEnumerationCount
             << ",\"runtime_unlimited\":true}" << endl;

        int index = -1;
        bool succeeded = true;
        ofstream sink;
        graphRepair::Result bound;
        const clock_t startedCpu = clock();
        const auto startedWall = chrono::steady_clock::now();
        if (method == "graph")
        {
            bound = graphRepair::calculate(molecule, false);
            index = bound.upperBound;
        }
        else
        {
            succeeded = improvedBnB(molecule, sink);
            index = lastCalculatedAssemblyIndex;
        }
        const double wallSeconds = chrono::duration<double>(
            chrono::steady_clock::now() - startedWall
        ).count();
        const double cpuSeconds = static_cast<double>(clock() - startedCpu)
            / static_cast<double>(CLOCKS_PER_SEC);
        const bool exactCompleted = method == "exact" && succeeded
            && !runtimeLimitReached && !enumerationLimitReached;
        cout << "{\"event\":\"finished\",\"method\":\"" << method
             << "\",\"algorithm_seconds\":" << wallSeconds
             << ",\"cpu_seconds\":" << cpuSeconds
             << ",\"assembly_index\":" << index
             << ",\"succeeded\":" << succeeded
             << ",\"exact_completed\":" << exactCompleted
             << ",\"runtime_limit_reached\":" << runtimeLimitReached
             << ",\"enumeration_limit_reached\":" << enumerationLimitReached
             << ",\"status\":\""
             << (method == "graph" ? "upper_bound" : exactCompleted ? "completed" : "incomplete")
             << "\"}" << endl;
        return succeeded ? 0 : 1;
    }
    catch (const exception &error)
    {
        cerr << "graph repair speed probe: " << error.what() << '\n';
        return 1;
    }
}
