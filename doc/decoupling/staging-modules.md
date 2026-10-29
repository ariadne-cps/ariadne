# Staging of System, Evolution and Analysis

The dependency graph of record is [structure.md](structure.md). This branch
installs six direct external dependencies in `ss/`: `symbolic`,
`integration`, `equation-solving`, `labelling`, `smt`, and `sat`.
`configuration` is pinned directly in `ss/cfg`.

The old local source modules `function/`, `geometry/`, `solving/`,
`symbolic/`, and `io/` were retired from `source/`. They are replaced by
the appropriate external repositories and their transitive dependencies;
there is deliberately no second local Symbolic implementation.

The only remaining Ariadne-owned source modules are:
- `dynamics/`, `hybrid/`: legacy holding areas for mixed responsibilities.
- `system/`, `evolution/`, `analysis/`: empty staging module targets.

The intended repository dependencies (dependent -> provider) are:

```text
analysis  -> evolution, smt, sat
evolution -> system, integration, labelling, equation-solving
system    -> symbolic
```

## Verification and migration blockers

- CMake dependency identity checks cover the major shared external gitlinks:
  Configuration, Symbolic, EquationSolving, Function, Cdcl, Bdd, Threading,
  and ConstraintSolving. These do not prove that legacy source headers are
  valid against the newly extracted repository APIs.
- `source/dynamics` and `source/hybrid` continue to refer to old
  `solving/` and `io/` headers. New upstream repositories place these APIs
  under `integration/`, `equation-solving/`, `smt/`, `sat/`, and the
  drawing/rendering chain. They cannot be assumed source-compatible.
- Existing Python aggregation, examples, tutorials, and packaging must be
  rewired for the repository structure; they are intentionally not configured
  as part of this staging pass.
- The legacy object libraries are retained but have NOT been proven to build.
  A passing CMake dependency consistency check is not a compile test.
- Before extraction, classify every source file by responsibility, update all
  includes, and compile public headers in isolation. Ensure `system` has no
  includes or target linkage to `evolution`, `analysis`, or `labelling`.
  Ensure `evolution` has no dependency on `analysis`.
- Restore actual compile/test, Python bindings, install/consumer and CI checks
  only after the relevant source/API migrations are complete.
