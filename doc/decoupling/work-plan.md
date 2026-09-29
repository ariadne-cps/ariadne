# Ariadne decoupling work plan

Status: active  
Branch: `decoupling`  
Baseline: `master` at `949d7e044ae65837fc02e6387701b10e7c15ddc6`  
Last updated: 2026-09-29

This is the living document for the decoupling activity. Update it in the same
commit as architectural changes that affect component boundaries, dependency
direction, packaging, or validation.

## Goal

Make the directories under `source/` independently buildable and testable with
declared, versionable dependencies. Separate repositories are an outcome of
verified boundaries; they are not the first mechanical step.

A component is ready for extraction when:

1. its public headers compile through an installed or exported target;
2. its CMake target declares every direct dependency with the correct
   `PUBLIC`, `PRIVATE`, or `INTERFACE` visibility;
3. its tests link to that target and to no aggregate Ariadne library;
4. a small external consumer can configure, compile, link, and run against the
   installed package;
5. dependency versions and configuration inputs are explicit;
6. no source include reaches into an undeclared sibling directory;
7. the component has no dependency cycle with a package intended to be
   released separately.

## Working principles

- Refactor and validate boundaries inside the monorepo before moving files.
- Preserve behaviour; structural commits should avoid unrelated API changes.
- Prefer moving misplaced integrations to a higher layer over introducing
  abstract interfaces without a concrete use case.
- Keep primitive mathematical types below algorithms that operate on them.
- Keep rendering adapters outside mathematical data types where practical.
- Treat header dependencies as part of the public contract, including template
  implementation headers.
- Add a focused validation before removing an include; transitive compilation
  success is not evidence of a sound public header.
- Record changes to the dependency graph and the rationale in this document.

## Target dependency shape

The precise package split remains a decision, but the intended direction is:

1. a small common foundation for utility-independent paradigms and logical
   contracts;
2. numeric types built on that foundation;
3. algebra built on numeric types;
4. primitive geometry built on numeric and algebra;
5. function and symbolic facilities organised around a shared low-level
   expression contract, without mutual implementation dependencies;
6. geometry algorithms, solvers, and rendering adapters above those primitives;
7. dynamics above the mathematical kernel;
8. hybrid above dynamics, with no dependency returning from dynamics or
   primitive geometry to hybrid.

This is a hypothesis to validate, not a frozen module list. In particular,
`geometry` currently mixes primitives, function-defined sets, paving
algorithms, solver calls, and drawing responsibilities. It is expected to split
before it becomes a standalone package.

## Milestones

### M0 — Reproducible baseline

- [x] Capture the direct include graph at the baseline commit.
- [x] Classify all direct relations as low, medium, or high.
- [x] Store the analysis, graph, summary, and include evidence in the repository.
- [ ] Add a repeatable repository script or build target that regenerates the
      dependency data and fails when an unreviewed edge appears.
- [ ] Record a clean baseline build and test result for each current aggregate:
      `ariadne-core`, `ariadne-kernel`, and `ariadne`.

Exit criterion: another contributor can regenerate the graph and reproduce the
baseline build without relying on this document's author.

### M1 — Explicit component targets

- [ ] Replace directory-wide include visibility with target include properties.
- [ ] Express every component-to-component dependency in CMake.
- [ ] Decide which headers are public, private, generated, or template
      implementation details.
- [ ] Generate `config.hpp` in the build tree and expose it through a target.
- [ ] Add public-header compilation checks for every component.
- [ ] Make each test directory link to its component target instead of an
      aggregate library.
- [ ] Register the currently omitted `tests/foundation` directory.

Exit criterion: removing an undeclared target dependency causes configuration
or compilation to fail locally and in CI.

### M2 — Remove small backward edges

- [ ] Verify and remove the apparently unused
      `dynamics/enclosure.cpp → hybrid/discrete_event.hpp` include.
- [ ] Verify and remove the apparently unused
      `geometry/list_set.hpp → hybrid/discrete_location.hpp` include.
- [ ] Isolate `algebra` drawing support for `Tensor` behind an adapter owned by
      the rendering layer.
- [ ] Move function-specific extensions from `algebra/graded.hpp` and
      `algebra/algebra_operations.tpl.hpp` to the function layer where feasible.
- [ ] Isolate the inclusion-integrator bridge to symbolic expression sets from
      the general solver interfaces.

Exit criterion: `hybrid` depends on `dynamics` and the mathematical kernel, but
neither `dynamics` nor lower-level geometry depends on `hybrid`.

### M3 — Stabilise the mathematical core

- [ ] Split generic logical/expression-template machinery from numeric and
      high-level symbolic specialisations.
- [ ] Resolve the `foundation ↔ numeric ↔ symbolic` cycles without moving the
      same cycle into a nominally lower package.
- [ ] Stabilise `algebra → function` as the primary direction by relocating the
      two function-specific algebra integrations.
- [ ] Identify the smallest primitive-geometry API needed by algebra and
      function (`Interval`, `Box`, declarations, and associated operations).
- [ ] Decide whether foundation and numeric initially ship as one repository
      with separate targets or as independently versioned packages.

Exit criterion: the foundation/numeric/algebra/primitive-geometry subgraph is
acyclic at package level and has standalone consumers.

### M4 — Separate mixed responsibilities

- [ ] Split geometry primitives from function-defined sets and paving
      algorithms.
- [ ] Move linear/nonlinear programming integrations out of primitive geometry.
- [ ] Separate graphics interfaces, mathematical drawing adapters, and
      Cairo/Gnuplot backends.
- [ ] Separate generic symbolic templates from conversions between
      `Expression`, `Formula`, functions, and sets.
- [ ] Re-evaluate component names and package boundaries after these moves.

Exit criterion: rendering and optimisation can be disabled without changing or
rebuilding the primitive mathematical packages.

### M5 — Repository extraction

- [ ] Select the first extraction candidate based on the verified graph.
- [ ] Provide install/export rules and a versioned CMake package.
- [ ] Add a standalone CI workflow and external-consumer test.
- [ ] Replace the monorepo source directory with a pinned dependency during a
      transition period.
- [ ] Document release compatibility and coordinated-change procedure.
- [ ] Repeat one component at a time; keep an integration build spanning all
      released packages.

Exit criterion: a component can be cloned, built, tested, installed, consumed,
and released without checking out the Ariadne monorepo.

## First work slice

The first implementation slice should be deliberately small and measurable:

1. add component target dependency declarations without moving files;
2. add public-header compilation checks for `foundation`, `numeric`, and
   `algebra`;
3. make their tests link to component-level targets;
4. remove the two suspected unused backward includes after compile and test
   verification;
5. regenerate the dependency graph and record the changed edge count.

Expected result: two low-cost cycle-closing edges disappear, while the build
starts enforcing the dependencies needed for the core work.

## Work log

| Date | Change | Validation | Graph impact | Status |
|---|---|---|---|---|
| 2026-09-29 | Replaced the in-tree `foundation` module with the standalone `ariadne-cps/foundation` submodule; removed local source/tests copies and wired `FOUNDATION_SRC` plus the `foundation` interface into Ariadne aggregates. | Verified `configuration@6b90c981939253a6740e93b726137ea3b6934ecb` is identical in Ariadne, Threading, Foundation and Utility; verified `utility@d28ec8bfa0f176f948756b00ac9c14e8eb3b2de2` is identical in Threading and Foundation. CI validation required after integration. | `foundation` is now an external repository boundary; local duplicate implementation removed. | In progress |
| 2026-09-29 | Made `foundation` source-level autonomous from `numeric` and `symbolic`: logical expression nodes are now owned by foundation; infinite `Sequence` conjunction/disjunction moved to `numeric/logical_sequence`; registered foundation tests and linked them only to `ariadne-foundation` + `utility`. | Static dependency assertions completed; build/test validation still required. | Expected removal of backward `numeric → foundation` and `symbolic → foundation`, breaking both mutual pairs while retaining intended `foundation → numeric`. | In progress |
| 2026-09-29 | Removed graphics and whole-tensor stream output from `algebra/Tensor`; moved Tensor drawing to the `io` layer through `tensor_drawable`, and adjusted the acoustic PDE example. | Static source review only; build/test validation still required. | Expected removal of backward `io → algebra`; the existing forward `algebra → io` remains. | In progress |
| 2026-09-29 | Removed the two low-volume backward includes from `geometry/list_set.hpp` and `dynamics/enclosure.cpp` into `hybrid`. | Static symbol inspection: neither consumer uses `DiscreteLocation`/`DiscreteEvent`; build and test validation still required. | Expected removal of `hybrid → geometry` and `hybrid → dynamics`; 55 → 53 direct edges pending regeneration. | In progress |
| 2026-09-29 | Created `decoupling` from `master`; recorded the baseline analysis and initial plan. | Branch SHA matched `master` before the documentation commit. | Baseline: 55 edges, 13 mutual pairs, one strongly connected component; 9 A / 35 M / 11 B. | Done |

## Decision log

| ID | Date | Decision | Rationale | Revisit when |
|---|---|---|---|---|
| D-001 | 2026-09-29 | Refactor boundaries in the monorepo before creating component repositories. | All ten directories are in one strongly connected component; immediate extraction would preserve cycles and add versioning overhead. | At least one component meets the extraction criteria. |
| D-002 | 2026-09-29 | Weight header propagation when prioritising coupling. | Header dependencies affect downstream consumers and templates, so raw include counts understate their cost. | A regeneration tool provides a better semantic metric. |
| D-003 | 2026-09-29 | Treat the proposed target dependency shape as provisional. | `geometry`, `symbolic`, and `io` mix responsibilities that must be separated before final package names are credible. | M3 and M4 produce tested boundaries. |

## Open questions

- Should foundation and numeric be released together initially while retaining
  separate targets?
- Which interval and box types belong to primitive geometry, and which aliases
  belong to function or set packages?
- Should drawing use free-function adapters, explicit renderer objects, or a
  small non-owning interface package?
- Where should the generic expression-template machinery live so that numeric,
  foundation, function, and symbolic do not form a cycle?
- Which component should be the first real repository extraction: hybrid as a
  high-level consumer, or a stable low-level core package?

## Risks

- Include removal may reveal accidental reliance on transitive headers.
- Template instantiation and link-time dependencies may not appear in the direct
  include graph.
- Splitting types across package boundaries can create ABI and release-lockstep
  constraints even after the include graph is acyclic.
- A large interface layer can disguise coupling instead of reducing it.
- Separate repositories can slow coordinated refactors unless an integration
  build continuously tests compatible revisions.

## Updating this document

For each decoupling change, update the work log with the validation command or
CI job, the dependency edges added or removed, and the relevant decision. Add a
new decision entry when a package boundary or dependency direction changes.
Keep completed milestones for historical context rather than deleting them.

