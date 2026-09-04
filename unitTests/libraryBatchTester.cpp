#include "assemblycpp.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// The implementation must not leak a generic global with this name from the
// installed static archive.
bool verbose = false;

namespace
{
class TemporaryDirectory
{
public:
    std::filesystem::path path;

    TemporaryDirectory()
    {
        const auto suffix = std::chrono::steady_clock::now()
                                .time_since_epoch()
                                .count();
        path = std::filesystem::temp_directory_path() /
               ("assemblycpp-library-test-" + std::to_string(suffix));
        std::filesystem::create_directory(path);
    }

    ~TemporaryDirectory()
    {
        std::error_code ignored;
        std::filesystem::remove_all(path, ignored);
    }
};

bool require(bool condition, const char *message)
{
    if (condition) return true;
    std::cerr << message << '\n';
    return false;
}
}

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        std::cerr <<
            "expected icosane, sucrose, butane, and native graph input paths\n";
        return 2;
    }

    const std::vector<std::string> inputs = {argv[1], argv[2]};
    const std::vector<assemblycpp::CalculationResult> batch =
        assemblycpp::calculateBatch(inputs);
    if (
        !require(batch.size() == 2, "batch result count mismatch") ||
        !require(batch[0].succeeded, "icosane batch calculation failed") ||
        !require(batch[0].assemblyIndex == 6, "icosane batch index mismatch") ||
        !require(batch[1].succeeded, "sucrose batch calculation failed") ||
        !require(batch[1].assemblyIndex == 8, "sucrose batch index mismatch")
    ) return 1;

    std::ifstream stream(argv[1]);
    const assemblycpp::CalculationResult streamed =
        assemblycpp::calculateMolfile(stream);
    if (
        !require(streamed.succeeded, "stream calculation failed") ||
        !require(streamed.assemblyIndex == 6, "stream calculation index mismatch")
    ) return 1;

    std::ifstream graphStream(argv[4]);
    const assemblycpp::CalculationResult streamedGraph =
        assemblycpp::calculateGraph(graphStream);
    const assemblycpp::CalculationResult graphFile =
        assemblycpp::calculate(argv[4]);
    if (
        !require(streamedGraph.succeeded, "graph stream calculation failed") ||
        !require(streamedGraph.input == "<stream>", "graph stream input mismatch") ||
        !require(
            streamedGraph.assemblyIndex == 5,
            "graph stream calculation index mismatch"
        ) ||
        !require(graphFile.succeeded, "graph file calculation failed") ||
        !require(
            streamedGraph.assemblyIndex == graphFile.assemblyIndex,
            "graph stream and file indices differ"
        )
    ) return 1;

    const std::string explicitHydrogenGraph =
        "explicit hydrogens\n"
        "6\n"
        "1 3 2 3 3 4 4 5 4 6\n"
        "H H C C H H\n"
        "1 1 1 1 1\n";
    std::istringstream filteredGraphStream(explicitHydrogenGraph);
    const assemblycpp::CalculationResult filteredGraph =
        assemblycpp::calculateGraph(filteredGraphStream);
    assemblycpp::CalculationOptions retainedHydrogenOptions;
    retainedHydrogenOptions.removeHydrogens = false;
    std::istringstream retainedGraphStream(explicitHydrogenGraph);
    const assemblycpp::CalculationResult retainedGraph =
        assemblycpp::calculateGraph(
            retainedGraphStream,
            retainedHydrogenOptions
        );
    if (
        !require(filteredGraph.succeeded, "filtered graph calculation failed") ||
        !require(filteredGraph.assemblyIndex == 0, "filtered graph index mismatch") ||
        !require(retainedGraph.succeeded, "retained graph calculation failed") ||
        !require(retainedGraph.assemblyIndex == 3, "retained graph index mismatch")
    ) return 1;

    std::istringstream invalidGraph(
        "invalid graph\n2\n1 3\nC C\n1\n"
    );
    const assemblycpp::CalculationResult rejectedGraph =
        assemblycpp::calculateGraph(invalidGraph);
    if (
        !require(!rejectedGraph.succeeded, "invalid graph stream succeeded") ||
        !require(
            rejectedGraph.error.find("outside the declared graph size") !=
                std::string::npos,
            "invalid graph stream omitted its parse error"
        )
    ) return 1;

    TemporaryDirectory temporaryDirectory;
    const std::filesystem::path copiedInput =
        temporaryDirectory.path / "icosane.mol";
    std::filesystem::copy_file(argv[1], copiedInput);
    const assemblycpp::CalculationResult noFileResult =
        assemblycpp::calculate(copiedInput.string());
    std::size_t fileCount = 0;
    for ([[maybe_unused]] const auto &entry :
         std::filesystem::directory_iterator(temporaryDirectory.path))
    {
        ++fileCount;
    }
    if (
        !require(noFileResult.succeeded, "no-file calculation failed") ||
        !require(noFileResult.assemblyIndex == 6, "no-file index mismatch") ||
        !require(fileCount == 1, "library calculation created an output file")
    ) return 1;

    assemblycpp::CalculationOptions limitedOptions;
    limitedOptions.runtimeTicks = 0;
    const assemblycpp::CalculationResult limited =
        assemblycpp::calculate(argv[3], limitedOptions);
    if (
        !require(limited.succeeded, "runtime-limited calculation failed") ||
        !require(limited.runtimeLimitReached, "runtime limit was not reported")
    ) return 1;

    const assemblycpp::CalculationResult afterLimit =
        assemblycpp::calculate(argv[3]);
    if (
        !require(afterLimit.succeeded, "post-limit calculation failed") ||
        !require(!afterLimit.runtimeLimitReached, "runtime stop leaked between calls") ||
        !require(afterLimit.assemblyIndex == 2, "post-limit index mismatch")
    ) return 1;

    assemblycpp::CalculationOptions invalidOptions;
    invalidOptions.enumerationLimit = 0;
    const assemblycpp::CalculationResult invalid =
        assemblycpp::calculate(argv[1], invalidOptions);
    if (
        !require(!invalid.succeeded, "invalid options unexpectedly succeeded") ||
        !require(!invalid.error.empty(), "invalid options omitted an error")
    ) return 1;

    const assemblycpp::CalculationResult missing =
        assemblycpp::calculate("assemblycpp-library-missing-input");
    if (
        !require(!missing.succeeded, "missing input unexpectedly succeeded") ||
        !require(!missing.error.empty(), "missing input omitted an error")
    ) return 1;

    return 0;
}
