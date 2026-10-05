# Installation

[Project overview](../README.md) · [Command line](cli.md) · [Build and tests](development.md)

Run source-tree examples from the repository root unless stated otherwise.

Install ParallelAssemblyCpp when you want a standalone command, reusable library,
and CMake package outside the build directory. With CMake 3.25 or newer, Ninja,
and a C++20 compiler available, run:

```bash
cmake --preset release
cmake --build --preset release
cmake --install build/release --prefix build/install
```

The first two commands can be skipped after completing the quick start.
`build/install` is a user-writable installation prefix, so administrator access
is not required. Verify the installed command with:

```bash
./build/install/bin/ParallelAssemblyCpp --help
```

The installation contains:

- The command-line tool in `<prefix>/bin`.
- The public header in `<prefix>/include/parallelassemblycpp`.
- The static library and CMake package files in the platform's library
  directory, typically `<prefix>/lib`.
- The README, license, and guides in `<prefix>/share/doc/parallelassemblycpp`.
  Development examples and fixtures require a source checkout.

Only the serial command and library are installed. The optional OpenMP, MPI,
hybrid, and telemetry executables are run from their build directory.

`--prefix` selects where the files are copied; it does not update `PATH`.
Replace `build/install` with another destination if needed, and add
`<prefix>/bin` to `PATH` to invoke `ParallelAssemblyCpp` from any directory.
Installing to a system location may require administrator privileges. On Windows,
the installed command is `build\install\bin\ParallelAssemblyCpp.exe`.

## Windows

The native Windows CI build uses Visual Studio 2022. In a Developer PowerShell
with CMake available, configure, build, and install the serial release:

```powershell
cmake -S . -B build/windows -G "Visual Studio 17 2022" -A x64 -DBUILD_TESTING=OFF
cmake --build build/windows --config Release
cmake --install build/windows --config Release --prefix build/install
.\build\install\bin\ParallelAssemblyCpp.exe --help
```

This generator does not require Ninja. The supplied Conda environment includes
Open MPI and is used for Linux development; it is not a native Windows setup.
Python is needed only if tests or benchmark scripts are enabled. The
`--memory-report=1` option is Linux-only.
