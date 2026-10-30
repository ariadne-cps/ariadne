# Staging of System, Evolution and Analysis

The dependency graph of record is [structure.md](structure.md). This branch
pins five direct external code dependencies in `ss/`: `integration`,
`equation-solving`, `labelled`, `smt`, and `sat`, plus the shared
`configuration` submodule in `ss/cfg`. `symbolic` is transitive, not a
direct Ariadne submodule.

The old local source modules `function/`, `geometry/`, `solving/`,
`symbolic/`, and `io/` were retired from `source/`. They are replaced by
the appropriate external repositories and their transitive dependencies;
there is deliberately no second local Symbolic implementation.

The Ariadne-owned source layers and their responsibilities are:

- `hybrid/`: fundamental hybrid representations, sets, paving, and graphics,
  analogous to the representations supplied by `labelled`.
- `system/`: system definitions, including discrete maps, vector fields,
  differential inclusions, and hybrid automata.
- `evolution/`: system evolution, including enclosures, evolvers, simulators,
  and orbits.
- `analysis/`: system analysis, including reachability and safety.
- `dynamics/`: legacy source location for some evolution and analysis code
  still awaiting physical relocation.

The current source-layer target dependencies (dependent -> provider) are:

```text
hybrid    -> labelled
system    -> hybrid
evolution -> system, integration, equation-solving
analysis  -> evolution, smt, sat
```

As a transitional exception, `DiscreteEvent` is still defined in
`source/system/discrete_event.hpp`, with a compatibility include in
`source/hybrid`; its intended ownership is the hybrid representation layer.
`DiscreteLocation` and `HybridSpace` already have canonical definitions in
`source/hybrid`. All three are domain primitives, not system algorithms.

## Verification and migration blockers

- CMake dependency identity checks cover the major shared external gitlinks:
  Configuration, Symbolic, EquationSolving, Function, Cdcl, Bdd, Threading,
  and ConstraintSolving. These do not prove that legacy source headers are
  valid against the newly extracted repository APIs.
- Legacy sources may still refer to old `solving/` and `io/` headers. The
  extracted repositories expose APIs through `integration/`,
  `equation-solving/`, `smt/`, `sat/`, and the drawing/rendering chain.
  Compile compatibility must be verified for each migrated source.
- Existing Python aggregation, examples, tutorials, and packaging must be
  rewired for the repository structure; they are intentionally not configured
  as part of this staging pass.
- Aggregate compilation and test CI were green for commit `25cfd9b` on
  Ubuntu, macOS, Windows, Coverage, and Sanitizers. This does not prove
  standalone public-header or module-boundary isolation.
- Before further extraction, classify sources by responsibility, update
  includes, and compile public headers in isolation. `system` may consume
  `hybrid` and its transitive `labelled` representation API, but must not
  depend on `evolution` or `analysis`; `evolution` must not depend on `analysis`.
- Restore actual compile/test, Python bindings, install/consumer and CI checks
  only after the relevant source/API migrations are complete.

## CI scope and configuration

Unix, Windows, Coverage, and Sanitizers run on `decoupling` as well as `main`.
They follow the extracted repositories by avoiding direct GMP/MPFR/Cairo
packages. Compilation and C++ tests of the existing `dynamics` and `hybrid`
modules remain mandatory. Python, installation and example/tutorial integration
have not yet been reconnected and remain explicit work items.

`configuration` is pinned to `cb75e7472b900025b7016d103c2e170197be1996`,
the latest verified `main` commit on 2026-10-29, shared by the dependencies.

## Local test ownership

Only `tests/dynamics` and `tests/hybrid` are retained as Ariadne-owned
C++ test suites. `tests/utility.hpp` is kept as their shared support header.
The former Function, Geometry, Solving, Symbolic and IO test directories and
obsolete `python/tests` were removed from this repository; those components
are tested in their owning external repositories. This does not exclude the
retained test executables from the CMake `tests` or CI targets.
