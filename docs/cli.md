# Command-line reference

[Project overview](../README.md) · [Command line](cli.md) · [Build and tests](development.md)

Run source-tree examples from the repository root unless stated otherwise.

```text
ParallelAssemblyCpp INPUT [OPTIONS]
ParallelAssemblyCpp [OPTIONS] -- INPUT
ParallelAssemblyCpp --help
```

By default, `INPUT` may be a V2000 MOL/SDF file or a ParallelAssemblyCpp native
graph file. With `--run-strings=1`, it is an exact text-file path instead. The
`.mol` and `.sdf` suffixes select MOL parsing case-insensitively. The exact path
is tried first; if it does not exist and has neither recognised suffix, a
lowercase `.mol` suffix is appended.
This allows `unitTests/alanine` to select `unitTests/alanine.mol`. An
`.sdf` input reads only the first record, which must be V2000; the rest of the
file is ignored. Native graph filenames must be supplied in full. Options may
appear before or after the input and use `--name=value` syntax. Boolean values
are `0` or `1`. Repeating an option, including through a legacy alias, is an
error.

Both graph formats are validated before the search starts. A bond is rejected
when it repeats an atom pair, joins an atom to itself, names an atom outside
the declared range, or has order zero, and an atom label is rejected when it is
not valid UTF-8. A rejected input produces a diagnostic on standard error and a
non-zero exit status, and writes no output files.

V2000 parsing uses the atom symbols and bond orders. Coordinates, charges,
isotopes, stereochemistry, and other property fields do not contribute to the
graph labels. The atom and bond blocks are read; trailing property blocks and
later SDF records are ignored.

A lone `--` ends option parsing, and every later argument is read as `INPUT`.
Use it to pass an input whose name begins with a dash, which is otherwise
rejected as an unknown option:

```bash
./build/release/ParallelAssemblyCpp -- -dashed-name.mol
```

Nothing after `--` is treated as an option, so a later `--help` or
`--name=value` is read as `INPUT` instead, and a second `INPUT` is an error.

## Native graph format

Supply all five lines, in this order:

1. A graph name.
2. The number of vertices.
3. Whitespace-separated, one-based endpoint pairs for all edges.
4. One whitespace-separated UTF-8 atom label per vertex.
5. One positive integer bond label per edge, in the same order as line 3.

For example, save the following as `chain.graph`:

```text
chain
3
1 2 2 3
C C C
1 1
```

Run `ParallelAssemblyCpp chain.graph`. The example describes three carbon
vertices joined by two single bonds. The vertex count may be zero; blank edge,
atom, or bond lines are still required when their lists are empty. Vertex counts and
bond labels must fit the parser's `short` representation (maximum 32767 on the
supported builds). Self-loops and repeated vertex pairs are rejected. Atom
labels cannot contain whitespace, and explicit `H` vertices are removed unless
`--remove-hydrogens=0`. Lines after the fifth are ignored.

## Options

| Option | Default | Purpose |
| --- | --- | --- |
| `-h`, `--help` | — | Show command help. |
| `--` | — | Stop parsing options; read every later argument as `INPUT`. |
| `--algorithm=<full\|re-pair>` | `full` | Select full exact search or a Re-Pair upper bound for graphs or strings. |
| `--runtime=<TICKS>` | Unlimited | Stop after the given `std::clock` budget. |
| `--enum-max=<COUNT>` | `50000000` | Limit retained connected masks in the initial graph DAG. |
| `--pathway=<0\|1>` | `1` | Write the recovered pathway. |
| `--run-strings=<0\|1>` | `0` | Treat `INPUT` as a file containing one string per line. |
| `--accept-palindromes=<0\|1>` | `0` | In string mode, identify a fragment with its reversal. |
| `--parallel=<auto\|on\|off>` | `off` | Select parallel search automatically, require it, or disable it. |
| `--threads=<auto\|N>` | `auto` | Set the OpenMP thread count per process; `N` must be positive. |
| `--remove-hydrogens=<0\|1>` | `1` | Remove explicit hydrogens from MOL/SDF and native graph inputs. |
| `--verbose=<0\|1>` | `0` | Print the parsed graph or each input string. |
| `--compensate-disjoint=<0\|1>` | `0` | For graphs, subtract one per processed component after the first. |
| `--upper-bound=graph-repair` | Disabled | Graph-only compatibility selector for `--algorithm=re-pair`; cannot be combined with `--algorithm`. |
| `--memory-report=<0\|1>` | `0` | Write Linux peak virtual memory to `memUsage`. |
| `--telemetry=<0\|1>` | `0` | Write search telemetry. |
| `--write-intermediate-mas=<0\|1>` | `0` | Write each improved index and its clock tick. |

`--enum-max` includes one-edge masks. Runtime and enumeration limits return the
best index found so far, which may not be the proven minimum. Run
`ParallelAssemblyCpp --help` for full details and accepted legacy option names.
`--telemetry` is available only in telemetry-enabled executables.
Graph mode rejects `--accept-palindromes=1`. String mode rejects explicit
`--enum-max` and `--remove-hydrogens` options, and rejects
`--compensate-disjoint=1`, because these options only apply to graphs.
Disabled boolean options such as `--compensate-disjoint=0` remain accepted.
`--threads` is only used when `--parallel` is `auto` or `on`.
`--memory-report=1` is rejected outside Linux, and fails if `memUsage` refers
to the input file, preserving the input.

`--runtime` is a cooperative budget in `std::clock` ticks, not a wall-clock
timeout. `CLOCKS_PER_SEC` is platform-dependent; it converts those ticks to the
clock's seconds. A zero budget returns the initial bound. In string mode the
budget applies independently to each line. Checks occur during the search, so
the budget is not a hard deadline for the entire process.

For molecular graphs, `--threads=auto` treats the OpenMP runtime default
(including `OMP_NUM_THREADS`) as an upper limit. Once the prepared root jobs
and DAG indicate enough work
for parallel search, it estimates one worker per 32,768 work units and rounds
the budget up to teams of 8, 16, 32, and so on. The eight-worker starting budget
leaves room for recursive work that this estimate can understate. The budget is
divided across launched MPI ranks, retaining at least one thread per rank.
Explicit thread counts apply to each process and are never reduced by this cap.

For molecular graphs in a parallel-enabled executable, `--parallel=auto`
prepares the root jobs and DAG, then uses their estimated search work to choose
parallel or serial execution. A serial fallback reports its reason.
`--parallel=on` forces parallel
search; automatic threads still use the workload cap, with at least two workers
in a single process. It fails if parallel execution cannot be honored, such as
when only one worker is available. `--parallel=off` always uses serial search.
Finite `--runtime` budgets and `--write-intermediate-mas=1` require serial
search, so `auto` reports a fallback and `on` reports an error. Pathway output
is supported with parallel optimization followed by bounded deterministic
reconstruction of a winning pathway.

## Outputs

For a MOL/SDF file, `INPUT` below excludes its recognised suffix.

| File | Contents |
| --- | --- |
| `INPUTOut` | Index or bound, `std::clock` ticks, and status lines for limited, interrupted, or heuristic calculations. |
| `INPUTPathway` | Graph pathway JSON when `--pathway=1`. |
| `INPUT_N_Pathway` | String pathway JSON for zero-based line `N` when `--pathway=1`. |
| `INPUTIntermediateMAs` | Elapsed ticks followed by the improved index, when `--write-intermediate-mas=1`. |
| `INPUTTelemetry.json` | Graph search counters when `--telemetry=1` in a telemetry-enabled executable. |
| `memUsage` | Linux `VmPeak` value when `--memory-report=1`. |

The pathway JSON is written as ASCII. Atom labels and string fragments outside
the ASCII range appear as `\uXXXX` escapes, so the file decodes the same way
whatever encoding the reader assumes.

Output paths are derived from the input; `memUsage` is in the working directory.
Enabled outputs overwrite files from previous runs. Outputs disabled on a later
run are not removed, so an older pathway or telemetry file can remain beside a
new `Out` file. Graph exact pathways, string exact pathways, and the two Re-Pair
certificates have different schemas; the Re-Pair formats are described below.

| Exit code | Meaning |
| --- | --- |
| `0` | Help shown or requested outputs completed, including limited searches and heuristic bounds. |
| `1` | Input, calculation, initialization, output, or execution-policy failure. |
| `2` | Command-line parser rejection, such as an unknown option or invalid value. |
| `130` | Cooperative user interruption after output handling succeeds. |

Successful exact runs do not add a separate completion status line. Check both
the exit code and any `status:` lines in `INPUTOut` before treating an index as
a proven minimum. An unavailable execution mode, such as `--parallel=on` in a
serial build, returns `1`. Some incompatible options are rejected by the parser
with `2`; others are rejected during execution-policy checks with `1`.
A failed output operation takes precedence over an interrupt
exit code, and partial output may exist after a failure. Invalid graph inputs
are rejected before any calculation output is opened.

## GraphRePair-inspired molecular upper bound

Select full exact search or the optional greedy bound with `--algorithm`.
Full search is the default; Re-Pair obtains a constructive assembly pathway
without enumerating every connected subgraph:

```bash
./build/release/ParallelAssemblyCpp molecule.mol --algorithm=full
./build/release/ParallelAssemblyCpp molecule.mol --algorithm=re-pair
```

Graph full mode runs the exact solver directly, without a Re-Pair prepass. The
graph Re-Pair calculation is available only when explicitly selected as
bound-only mode. The same distinction applies to the library's
`graphRepairUpperBound` option.

The earlier `--upper-bound=graph-repair` option is still supported. Use one
selector per command; combining it with `--algorithm` is an error.

The algorithm repeatedly combines adjacent fragments whose labelled graph occurs
more than once, sharing the construction across edge-disjoint occurrences.
Its bound counts the binary joins needed to build the reusable fragments and
assemble the remaining graph. Atom labels, bond orders, and all attachment
vertices are retained. Explicit hydrogen removal follows the usual
`--remove-hydrogens` setting.

`INPUTOut` reports the heuristic bound, the trivial bound, and the status
`heuristic upper bound (minimum not proven)`. `INPUTPathway` contains a
`graph-repair-assembly-v1` JSON certificate with reusable construction rules and
the residual assembly; this format differs from the exact solver's pathway
format. `--pathway=0` skips that file. A smaller bound means a shorter known
construction, and does not establish that the minimum has been found.

The default remains exact search. The heuristic is serial: `--parallel=auto`
falls back and `--parallel=on` is rejected. Explicit `--runtime` and
`--enum-max`, and enabled telemetry or intermediate-index output are
unavailable in this mode. String inputs use the string counterpart described
below; the legacy `--upper-bound=graph-repair` selector remains graph-only.
For disconnected graphs, `--compensate-disjoint=1`
subtracts the joins between components that contain bonds. Isolated atoms cost
no joins, and an input with no bonds has bound zero.

## String assembly

Pass `--run-strings=1` to select the string algorithm. In this mode `INPUT` is
opened exactly as supplied and every line is processed as a separate string:

```bash
./build/release/ParallelAssemblyCpp strings.txt --run-strings=1
```

String files must contain valid UTF-8. Each Unicode code point is one symbol;
no Unicode normalization is applied. LF and CRLF line endings are accepted,
empty lines are separate strings, and the last line need not end in a newline.
Pathway positions and lengths count code points rather than UTF-8 bytes.
The empty string has assembly index `-1`; a one-symbol string has index `0`.

Exact string search starts with a Re-Pair construction and retains its pathway
as the initial incumbent. LZ-style bounds then prune branches that cannot
improve it, including ties; the retained witness remains valid when the seed
is already optimal. Search can improve a suboptimal seed and proves the
minimum on completion. The prepass shares the search's runtime and cancellation
budget; a zero budget or immediate cancellation retains the trivial bound.

Results are written to `strings.txtOut`. With pathway output enabled, the
zero-based line number is included in each pathway name, such as
`strings.txt_0_Pathway`. `--accept-palindromes=1` treats a fragment and its
reversal as equivalent. OpenMP, MPI, and hybrid executables can distribute
the search branches within each string using the same `--parallel` and
`--threads` options as molecular search:

```bash
./build/parallel/ParallelAssemblyCppOMP strings.txt \
  --run-strings=1 --parallel=on --threads=4
```

Lines are still read and written in input order. Each worker has private search
caches seeded from the same root enumeration; OpenMP workers share the best
index, and MPI ranks search disjoint sets of root branches. Completed searches
reconstruct the pathway in serial order, so indices and pathway JSON match
serial execution. Use `--pathway=0` to skip reconstruction. Parallel overhead
can outweigh the benefit for short strings or searches with few branches;
OpenMP teams are capped by the number of local root branches.
With compatible options and multiple workers, `--parallel=auto` uses parallel
string search without applying the molecular DAG workload threshold.

Telemetry and intermediate-index output remain unavailable in string mode.
A finite `--runtime` budget applies separately to each line and requires serial
search: `--parallel=auto` falls back, while `--parallel=on` reports an error.
The serial executable rejects `--parallel=on`.

### String Re-Pair upper bound

Use the same algorithm selector to build a reusable binary string grammar:

```bash
./build/release/ParallelAssemblyCpp strings.txt --run-strings=1 --algorithm=re-pair
./build/release/ParallelAssemblyCpp strings.txt --run-strings=1 \
  --algorithm=re-pair --accept-palindromes=1
```

Adjacent fragments are grouped by their full expanded string, with reversal
equivalence when requested. Each round chooses the group with the greatest
saving across nonoverlapping occurrences; ties choose the earliest occurrence.
Creating a binary production costs one join and using an existing production
costs zero. The final bound is the number of productions plus the number of
residual fragments minus one. For example, `abababab` has bound 3, obtained by
building `ab`, building `abab`, then joining two copies. Empty and one-symbol
strings retain bounds -1 and 0.

`INPUTOut` explicitly reports an assembly upper bound, the status
`heuristic upper bound (minimum not proven)`, the trivial bound, rule count,
remaining fragments, and elapsed ticks for each line. With `--pathway=1`, each
`INPUT_N_Pathway` contains a `string-repair-assembly-v1` certificate. It records
the original string, Unicode scalar terminals, topologically ordered binary
rules with child orientations, and ordered residual occurrences with scalar
offsets, lengths, and orientations. Expanding those occurrences reconstructs
the input and independently verifies the construction cost. This certificate
format differs from the exact string pathway format. `--pathway=0` skips it.

String Re-Pair is serial, including in OpenMP and MPI executables.
`--parallel=auto` reports a serial fallback; `--parallel=on`, explicit runtime
or enumeration limits, telemetry, and intermediate-index output are rejected.
UTF-8 validation, line handling, reversal semantics, and output filenames match
full string mode. `--algorithm=full` remains the default and uses the construction
to seed exact search. A Re-Pair bound alone may exceed the minimum.

The implementation uses collision-free substring ranks and updates adjacent
occurrences locally, with deterministic frequency selection. Time and space
are O(n log n) for n Unicode scalars, including substring-rank preprocessing;
expanded fragments are not repeatedly copied or compared.

The implementation is adapted from the standard-library string search in
[AssemblyCPP Public](https://gitlab.com/croningroup/public/assemblycpp-public)
at commit `2a87948`, authored by Stuart Marshall from work by Ian Seet and
Leroy Cronin. It ports the interval enumeration, Lempel-Ziv lower bound, and
recursive branch-and-bound search without importing the upstream Boost/VF2
graph implementation or vendored dependencies. The legacy spellings
`-runStrings`, `-acceptPalindromes`, and `-palindrome` remain accepted.
