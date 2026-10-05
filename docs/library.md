# C++ library

[Project overview](../README.md) · [Command line](cli.md) · [Build and tests](development.md)

Installed packages export `ParallelAssemblyCpp::Library`. After
[installing the library](installation.md), create a separate consumer directory
with this `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.25)
project(AssemblyExample LANGUAGES CXX)
find_package(ParallelAssemblyCpp 0.1.0 CONFIG REQUIRED)
add_executable(my_program main.cpp)
target_link_libraries(my_program PRIVATE ParallelAssemblyCpp::Library)
```

Save this as `main.cpp` in the same directory:

```cpp
#include <parallelassemblycpp.h>

#include <iostream>

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: my_program INPUT\n";
        return 2;
    }
    const auto result = parallelassemblycpp::calculate(argv[1]);
    if (!result)
    {
        std::cerr << result.error << '\n';
        return 1;
    }
    std::cout << result.assemblyIndex << '\n';
    if (result.runtimeLimitReached || result.enumerationLimitReached || result.upperBoundOnly)
        std::cout << "Upper bound; minimum not proven\n";
}
```

Configure and build from the consumer directory, substituting the actual
installation prefix and input path:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/absolute/path/to/parallelassemblycpp/build/install
cmake --build build
./build/my_program /absolute/path/to/parallelassemblycpp/unitTests/alanine.mol
```

The imported target supplies the header path and C++20 requirement. With a
multi-configuration generator, use `--config Release` when building and run the
executable from `build/Release` (with `.exe` on Windows).

`calculateMolfile` accepts a V2000 molfile stream, while `calculateGraph`
accepts a ParallelAssemblyCpp native graph stream. `calculateBatch` processes
several inputs sequentially without process startup between items. Library calls
are serial even when parallel CLI targets are enabled. They do not create output
files or return pathway JSON. Search state is process-global, so the API is
reusable but not thread-safe; use separate processes for concurrent work.

A successful `CalculationResult` means a calculation produced an index; it does
not by itself establish minimality. The result is a proven minimum only when
`runtimeLimitReached`, `enumerationLimitReached`, and `upperBoundOnly` are all
false. For graph calls, `clockTicks` excludes input parsing. Batch results
preserve input order, and a failed or limited item does not stop later items.

Set `CalculationOptions::graphRepairUpperBound = true` to select the experimental
GraphRePair-inspired bound through `calculate`, `calculateMolfile`,
`calculateGraph`, or `calculateBatch`. A successful result
then has `CalculationResult::upperBoundOnly == true`, and `assemblyIndex` holds
the upper bound. The runtime budget must remain unlimited; the enumeration limit
must still be positive but does not limit this mode. Exact calls retain the
default `upperBoundOnly == false` (runtime or enumeration limits can still
prevent an exact proof).

`calculateString` and `calculateStringBatch` accept literal UTF-8 strings without
creating files. They use `StringCalculationOptions`: `acceptReversed` enables
reversal equivalence, `runtimeTicks` limits exact search, and `rePairUpperBound`
selects string Re-Pair with an unlimited runtime budget. Results use the same
`CalculationResult` status fields, with `input == "<string>"`; invalid UTF-8 is
reported in `error`, and batches continue after failed items. Embedded NULs and
newlines are ordinary symbols in these literal-string entry points. String
`clockTicks` includes decoding.

```cpp
parallelassemblycpp::StringCalculationOptions options;
options.rePairUpperBound = true;
const auto bound = parallelassemblycpp::calculateString("abababab", options);
// bound.assemblyIndex == 3; bound.upperBoundOnly == true
```
