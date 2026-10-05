# Algorithm and provenance

[Project overview](../README.md) · [Command line](cli.md)

## Molecular assembly

The molecular search implements the algorithm described by Ian Seet, Keith Y.
Patarroyo, Gage Siebert, Sara I. Walker, and Leroy Cronin in
[*Rapid Exploration of Assembly Chemical Space of Molecular Graphs*](https://arxiv.org/abs/2410.09100).

ParallelAssemblyCpp treats a molecule as a labelled graph: atoms are vertices and
bonds are edges. Its assembly index is the smallest number of joining steps
needed to build that graph when a fragment that has already been made can be
reused.

The program parses a V2000 MOL/SDF or native graph file, removes explicit
hydrogens by default, and enumerates connected fragments of the molecular
graph. It canonicalises those fragments so that structurally equivalent copies
can be recognised even when they use different atom or bond indices. A compact
directed acyclic graph (DAG) records the fragment relationships for reuse in
later search passes.

A branch-and-bound search then explores ways to reuse matching, disjoint
fragments. Lower bounds and a cache of previously visited canonical assembly
states eliminate branches that cannot improve the best result. The final
output contains the lowest assembly index found and, unless disabled, a JSON
description of a corresponding assembly pathway. If a runtime or enumeration
limit is reached, the reported index is the best result found so far and may
not be the proven minimum.

The [Re-Pair modes](cli.md#graphrepair-inspired-molecular-upper-bound) build
constructive upper bounds without proving optimality. Exact molecular search
does not use a Re-Pair prepass. Exact string search does use a Re-Pair
construction as its initial incumbent; see [string assembly](cli.md#string-assembly).
GPU backends are not available build targets.

## Compared with the original AssemblyCpp v5

This repository and the
[original AssemblyCpp v5](https://github.com/croningp/assemblycpp-v5) implement
the same molecular-assembly calculation. Both parse a labelled molecular graph,
enumerate connected fragments, identify structurally equivalent copies, search
possible fragment reuses with branch-and-bound, and recover a pathway for the
lowest assembly index found. This is therefore a re-engineering of the v5
implementation, not a different definition of the assembly index.

The comparison below is against the original repository's `main` branch at
[commit `f9209034`](https://github.com/croningp/assemblycpp-v5/commit/f9209034b0851d03282322bb6be697beaf030dda):

- **Graph representation and matching.** The original uses fixed 512-bit edge
  masks and Boost's VF2 implementation for cyclic graph isomorphism. This
  version uses compact, dynamically sized edge masks and in-project tree and
  cyclic canonicalisation, removing the vendored Boost dependency and the
  fixed 512-edge mask limit.
- **Search implementation.** This version adds compact DAG storage with
  immutable retained node masks, frontier-driven enumeration, reusable
  canonical fragment identities, residual and transposition caches, tighter
  bounds, and allocation reuse.
  These changes reduce repeated graph and search work without changing the
  quantity being calculated.
- **Parallel execution.** The original solver is serial. Serial search remains
  the default here, while optional OpenMP, MPI, and hybrid builds can distribute
  independent search branches between workers and then deterministically
  reconstruct a winning pathway. This can reduce runtime for larger graphs,
  where the search exposes enough work to outweigh parallel coordination
  overhead; small graphs may see little or no speed-up.
- **Interfaces.** The original provides a C++17 `assembly` command with
  file-based results. This project provides a C++20 `ParallelAssemblyCpp` command
  plus an installable `ParallelAssemblyCpp::Library` with stream, file, and batch
  APIs that return results directly. Its command-line handling also validates
  options and reports interrupted or limited searches explicitly.
- **Project tooling.** This version expands the build and verification support
  with CMake presets, package installation and export, CI, focused and full
  regression suites, pathway and parallel-parity tests, telemetry, maintained
  benchmarks, and optional LTO and PGO builds.

When allowed to finish, both implementations target the same minimum assembly
index. A recovered pathway can differ when more than one optimal pathway
exists, and a runtime or enumeration limit can make this version return a
best-so-far result instead of a proven minimum.
