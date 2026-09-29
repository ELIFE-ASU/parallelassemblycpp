#include "parallelassemblycpp.h"

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
               ("parallelassemblycpp-library-test-" + std::to_string(suffix));
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
    const std::vector<parallelassemblycpp::CalculationResult> batch =
        parallelassemblycpp::calculateBatch(inputs);
    if (
        !require(batch.size() == 2, "batch result count mismatch") ||
        !require(batch[0].succeeded, "icosane batch calculation failed") ||
        !require(batch[0].assemblyIndex == 6, "icosane batch index mismatch") ||
        !require(batch[1].succeeded, "sucrose batch calculation failed") ||
        !require(batch[1].assemblyIndex == 8, "sucrose batch index mismatch")
    ) return 1;

    std::ifstream stream(argv[1]);
    const parallelassemblycpp::CalculationResult streamed =
        parallelassemblycpp::calculateMolfile(stream);
    if (
        !require(streamed.succeeded, "stream calculation failed") ||
        !require(streamed.assemblyIndex == 6, "stream calculation index mismatch")
    ) return 1;

    std::ifstream graphStream(argv[4]);
    const parallelassemblycpp::CalculationResult streamedGraph =
        parallelassemblycpp::calculateGraph(graphStream);
    const parallelassemblycpp::CalculationResult graphFile =
        parallelassemblycpp::calculate(argv[4]);
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
    const parallelassemblycpp::CalculationResult filteredGraph =
        parallelassemblycpp::calculateGraph(filteredGraphStream);
    parallelassemblycpp::CalculationOptions retainedHydrogenOptions;
    retainedHydrogenOptions.removeHydrogens = false;
    std::istringstream retainedGraphStream(explicitHydrogenGraph);
    const parallelassemblycpp::CalculationResult retainedGraph =
        parallelassemblycpp::calculateGraph(
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
    const parallelassemblycpp::CalculationResult rejectedGraph =
        parallelassemblycpp::calculateGraph(invalidGraph);
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
    const parallelassemblycpp::CalculationResult noFileResult =
        parallelassemblycpp::calculate(copiedInput.string());
    parallelassemblycpp::CalculationOptions boundOptions;
    boundOptions.graphRepairUpperBound = true;
    const parallelassemblycpp::CalculationResult boundFileResult =
        parallelassemblycpp::calculate(copiedInput.string(), boundOptions);
    std::size_t fileCount = 0;
    for ([[maybe_unused]] const auto &entry :
         std::filesystem::directory_iterator(temporaryDirectory.path))
    {
        ++fileCount;
    }
    if (
        !require(noFileResult.succeeded, "no-file calculation failed") ||
        !require(noFileResult.assemblyIndex == 6, "no-file index mismatch") ||
        !require(!noFileResult.upperBoundOnly, "exact result marked as heuristic") ||
        !require(boundFileResult.succeeded, "bound file calculation failed") ||
        !require(boundFileResult.upperBoundOnly, "bound result omitted heuristic status") ||
        !require(
            boundFileResult.assemblyIndex >= 6 && boundFileResult.assemblyIndex < 18,
            "icosane bound is invalid or does not improve its trivial bound"
        ) ||
        !require(fileCount == 1, "library calculation created an output file")
    ) return 1;

    std::ifstream boundMolfileStream(argv[1]);
    const auto boundMolfile = parallelassemblycpp::calculateMolfile(
        boundMolfileStream, boundOptions
    );
    const auto boundBatch = parallelassemblycpp::calculateBatch(inputs, boundOptions);
    if (
        !require(boundMolfile.succeeded, "bound molfile stream failed") ||
        !require(boundMolfile.upperBoundOnly, "bound molfile omitted heuristic status") ||
        !require(
            boundMolfile.assemblyIndex == boundFileResult.assemblyIndex,
            "bound molfile stream and file disagree"
        ) ||
        !require(boundBatch.size() == 2, "bound batch result count mismatch") ||
        !require(
            boundBatch[0].succeeded && boundBatch[1].succeeded &&
            boundBatch[0].upperBoundOnly && boundBatch[1].upperBoundOnly,
            "bound batch failed or omitted heuristic status"
        ) ||
        !require(
            boundBatch[0].assemblyIndex == boundFileResult.assemblyIndex &&
            boundBatch[1].assemblyIndex >= batch[1].assemblyIndex,
            "bound batch produced an invalid bound"
        )
    ) return 1;

    std::istringstream boundHydrogenStream(explicitHydrogenGraph);
    const auto filteredBound = parallelassemblycpp::calculateGraph(
        boundHydrogenStream, boundOptions
    );
    boundOptions.removeHydrogens = false;
    std::istringstream boundRetainedStream(explicitHydrogenGraph);
    const auto retainedBound = parallelassemblycpp::calculateGraph(
        boundRetainedStream, boundOptions
    );
    if (
        !require(
            filteredBound.succeeded && filteredBound.upperBoundOnly &&
                filteredBound.assemblyIndex == 0,
            "bound mode failed to remove explicit hydrogen"
        ) ||
        !require(
            retainedBound.succeeded && retainedBound.upperBoundOnly &&
                retainedBound.assemblyIndex >= 3 && retainedBound.assemblyIndex <= 4,
            "bound mode failed to retain explicit hydrogen"
        )
    ) return 1;

    const std::string disconnectedGraph =
        "disconnected and isolated\n5\n1 2 3 4\nC C C C N\n1 1\n";
    std::istringstream disconnectedDefaultStream(disconnectedGraph);
    const auto defaultBound = parallelassemblycpp::calculateGraph(
        disconnectedDefaultStream, boundOptions
    );
    boundOptions.compensateDisjoint = true;
    std::istringstream disconnectedCompensatedStream(disconnectedGraph);
    const auto compensatedBound = parallelassemblycpp::calculateGraph(
        disconnectedCompensatedStream, boundOptions
    );
    std::istringstream isolatedStream("isolated\n1\n\nN\n\n");
    const auto isolatedBound = parallelassemblycpp::calculateGraph(
        isolatedStream, boundOptions
    );
    if (
        !require(
            defaultBound.succeeded && defaultBound.assemblyIndex == 1,
            "default disconnected bound mismatch"
        ) ||
        !require(
            compensatedBound.succeeded && compensatedBound.assemblyIndex == 0,
            "compensated disconnected bound counts isolated atoms"
        ) ||
        !require(
            isolatedBound.succeeded && isolatedBound.upperBoundOnly &&
                isolatedBound.assemblyIndex == 0,
            "empty-bond bound must be zero"
        )
    ) return 1;

    boundOptions.runtimeTicks = 0;
    const auto invalidBoundBudget = parallelassemblycpp::calculate(argv[3], boundOptions);
    const auto exactAfterBound = parallelassemblycpp::calculate(argv[3]);
    if (
        !require(
            !invalidBoundBudget.succeeded && !invalidBoundBudget.error.empty(),
            "bound mode accepted an unsupported runtime budget"
        ) ||
        !require(
            exactAfterBound.succeeded && !exactAfterBound.upperBoundOnly &&
                exactAfterBound.assemblyIndex == 2,
            "bound options leaked into subsequent exact calculation"
        )
    ) return 1;

    parallelassemblycpp::CalculationOptions limitedOptions;
    limitedOptions.runtimeTicks = 0;
    const parallelassemblycpp::CalculationResult limited =
        parallelassemblycpp::calculate(argv[3], limitedOptions);
    if (
        !require(limited.succeeded, "runtime-limited calculation failed") ||
        !require(limited.runtimeLimitReached, "runtime limit was not reported")
    ) return 1;

    const parallelassemblycpp::CalculationResult afterLimit =
        parallelassemblycpp::calculate(argv[3]);
    if (
        !require(afterLimit.succeeded, "post-limit calculation failed") ||
        !require(!afterLimit.runtimeLimitReached, "runtime stop leaked between calls") ||
        !require(afterLimit.assemblyIndex == 2, "post-limit index mismatch")
    ) return 1;

    // Upper-bound mode must not seed a subsequent exact calculation.
    const std::string eightBondChain =
        "eight-bond chain\n9\n"
        "1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9\n"
        "C C C C C C C C C\n1 1 1 1 1 1 1 1\n";
    parallelassemblycpp::CalculationOptions chainBoundOptions;
    chainBoundOptions.graphRepairUpperBound = true;
    std::istringstream fastChainStream(eightBondChain);
    const auto fastChain = parallelassemblycpp::calculateGraph(
        fastChainStream, chainBoundOptions
    );
    parallelassemblycpp::CalculationOptions chainLimitedOptions;
    chainLimitedOptions.runtimeTicks = 0;
    std::istringstream limitedChainStream(eightBondChain);
    const auto limitedChain = parallelassemblycpp::calculateGraph(
        limitedChainStream, chainLimitedOptions
    );
    std::istringstream exactChainStream(eightBondChain);
    const auto exactChain = parallelassemblycpp::calculateGraph(exactChainStream);
    if (
        !require(
            fastChain.succeeded && fastChain.upperBoundOnly &&
                fastChain.assemblyIndex == 3,
            "chain upper-bound calculation failed"
        ) ||
        !require(
            limitedChain.succeeded && !limitedChain.upperBoundOnly &&
                limitedChain.runtimeLimitReached && limitedChain.assemblyIndex == 7,
            "limited exact calculation inherited a Re-Pair seed"
        ) ||
        !require(
            exactChain.succeeded && !exactChain.upperBoundOnly &&
                !exactChain.runtimeLimitReached && exactChain.assemblyIndex == 3,
            "chain exact calculation inherited upper-bound or runtime options"
        )
    ) return 1;

    parallelassemblycpp::CalculationOptions invalidOptions;
    invalidOptions.enumerationLimit = 0;
    const parallelassemblycpp::CalculationResult invalid =
        parallelassemblycpp::calculate(argv[1], invalidOptions);
    if (
        !require(!invalid.succeeded, "invalid options unexpectedly succeeded") ||
        !require(!invalid.error.empty(), "invalid options omitted an error")
    ) return 1;

    const parallelassemblycpp::CalculationResult missing =
        parallelassemblycpp::calculate("parallelassemblycpp-library-missing-input");
    if (
        !require(!missing.succeeded, "missing input unexpectedly succeeded") ||
        !require(!missing.error.empty(), "missing input omitted an error")
    ) return 1;

    return 0;
}
