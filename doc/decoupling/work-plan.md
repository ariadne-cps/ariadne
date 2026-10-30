# Ariadne decoupling work plan

Status: active  
Base state: low-level chain through Algebra extracted and green  
Last updated: 2026-10-30

## Established chain

```text
utility -> foundation -> numeric -> interval -> algebra -> function
```

Utility, Foundation, Numeric, Interval and Algebra are standalone. The completed
Foundation/Algebra extraction history is recorded in Git history and
`current-state.md`; completed migration checklists are intentionally omitted
from this active plan.

## Target repository structure

The current source-layer graph is described in `structure.md`. In Ariadne the
responsibilities and dependency direction are:

```text
hybrid    -> labelled
system    -> hybrid
evolution -> system, integration, equation-solving
analysis  -> evolution, smt, sat
ariadne   -> analysis
```

`hybrid`, `system`, `evolution`, and `analysis` are Ariadne-owned source modules.
`labelled`, `integration`, `equation-solving`, `smt`, and `sat` are external
repositories. Dependencies such as `symbolic` are obtained transitively.

The ownership rule is semantic:

- `hybrid` owns fundamental hybrid-domain representations, including
  `DiscreteLocation`, `HybridSpace`, hybrid sets, paving,
  and graphics, just as `labelled` owns labelled geometric representations;
- `system` defines systems: discrete maps, vector fields, differential
  inclusions, discrete events, and hybrid automata;
- `evolution` owns how such systems evolve: enclosures, evolvers, simulators,
  and orbit types, including their hybrid variants;
- `analysis` owns how systems and their evolution are analysed: global
  storage, reachability, safety certificates, and related algorithms.

`system` consumes `hybrid` but does not own hybrid primitives. `evolution` and
`analysis` must not leak their algorithms into `system` or `hybrid`.
Transitive access to `labelled` via `hybrid` is intentional.

## F1 - Remove accidental Function infrastructure dependencies

- [ ] Remove the unused `logging/logging.hpp` include from
      `function/calculus_base.hpp`.
- [ ] Verify that Function has no semantic Logging or Threading dependency.
- [ ] Make Function target dependencies explicit rather than relying on
      Ariadne-wide linkage.

## F2 - Remove Function -> high-level Geometry

- [ ] Identify the minimal Box/domain value types required by Function.
- [ ] Separate those primitives from Point/SetInterface and higher-level
      Geometry responsibilities.
- [ ] Move measurable-function and multifunction set integrations above the
      Function boundary where appropriate.
- [ ] Ensure Function public headers no longer depend on the high-level
      `geometry/` package.

## F3 - Remove Function -> high-level Symbolic

- [ ] Separate generic symbolic machinery genuinely required by Function from
      high-level Expression/Space/Variable integration.
- [ ] Move concrete Expression/Function conversion bridges to the higher layer.
- [ ] Ensure Function no longer depends on high-level Symbolic implementation.

## F4 - Make Function standalone

- [ ] Define the standalone Function target above Algebra using only intended
      lower dependencies.
- [ ] Register Function public headers and use the dependency-bundle
      installation model.
- [ ] Give Function a direct Configuration dependency and verify matching
      Configuration gitlinks with its direct dependency.
- [ ] Add public-header isolation checks.
- [ ] Add focused C++ tests.
- [ ] Layer Python bindings on the aggregate `pyariadne-algebra` interface
      only.
- [ ] Add standalone installation and external-consumer/tutorial checks.
- [ ] Get Unix, Windows, Coverage and Python CI green.

## F5 - Integrate standalone Function into Ariadne

- [ ] Add the standalone `ariadne-cps/function` repository as Ariadne's direct
      lower Function dependency.
- [ ] Remove Ariadne's local Function implementation and Function-owned tests or
      bindings that move with the repository.
- [ ] Make Ariadne consume the aggregate Function C++ and Python interfaces.
- [ ] Propagate compatible dependency pins.
- [ ] Get the full Ariadne integration CI green.

## H1 - Establish the Labelled boundary

`labelled` is a normal repository in the dependency DAG, not a transverse
layer. Its responsibility is the generic labelled forms of concepts owned by
`kernel`.

- [ ] Move generic labelled forms of Kernel concepts into `labelled`.
      Candidates include `LabelledPoint`, `LabelledSet`, `LabelledBox`,
      `LabelledInterval`, and other labelled wrappers of Geometry/Function
      concepts.
- [ ] Classify `LabelledGrid`, `LabelledGridTreePaving`,
      `LabelledMapping`, labelled curves, and graphics-related labelled types
      according to the same rule: they belong in `labelled` only when the
      underlying unlabelled concept is owned by `kernel`.
- [ ] Keep `RealSpace`, variables, assignments, expressions, predicates,
      `VariablesPoint`, `VariablesBox`, and `RealVariablesBox` in
      `symbolic`. A named symbolic space is not by itself a labelled wrapper
      around a Kernel object.
- [ ] Split the current `symbolic/expression_set.hpp` responsibilities so that
      generic symbolic variable/domain machinery remains in `symbolic` while
      `LabelledSet`/`LabelledBox`-style wrappers move to `labelled`.
- [ ] Ensure `labelled` depends on `symbolic` and obtains `kernel`
      transitively through it.

The naming rule is deliberately not "every class beginning with
`Labelled` belongs to `labelled`". A thin labelled adaptation of a concept
introduced by a higher repository stays with that concept.

## H2 - Extract System definitions above Hybrid

`system` defines what a system is, not its evolution or analysis. The
fundamental hybrid representation layer stays in `source/hybrid`.

- [ ] Keep `DiscreteLocation`, `HybridSpace`, and other fundamental hybrid
      representations in `source/hybrid`, alongside hybrid sets, paving,
      and graphics. Do not relocate them to `system`.
- [x] Keep `DiscreteEvent` in `source/system`, because events belong to a
      system's transitions and behaviour rather than the hybrid state space.
- [ ] Consolidate system-definition classes in `source/system`:
      `IteratedMap`, `VectorField`, `DifferentialInclusion`,
      `HybridAutomaton`, `HybridAutomatonInterface`, and directly related
      specification types.
- [ ] Remove `EvolverType`, `EnclosureType`, and `StorageType` aliases from
      `IteratedMap`, `VectorField`, `DifferentialInclusion`, and
      `HybridAutomatonInterface`.
- [ ] Keep evolution, enclosure, and storage implementations out of System
      public headers. System may use lower-layer hybrid representations.
- [ ] Express system state variables using symbolic concepts (`RealSpace`,
      variables, assignments, and predicates); preserve
      `DifferentialInclusion` input domains expressed through
      `RealVariablesBox`.
- [ ] Keep analysis-specific space traits and set-approximation choices in
      `analysis`, without removing the hybrid sets and graphics owned by
      `hybrid`.
- [ ] Remove actual unwanted dependencies on `evolution` or `analysis`
      from System headers and sources; keep only the ordinary module build.
- [ ] Validate the resulting dependencies through the existing C++ build and
      tests, without adding special CMake header-isolation targets.

## H3 - Extract Evolution

`evolution` owns local representations of evolving sets and the algorithms
that advance a system.

- [ ] Move `Enclosure`, evolver interfaces and implementations, simulators,
      and orbit machinery into `evolution`.
- [ ] Move the corresponding hybrid evolution types into the same repository,
      including `HybridEnclosure`, hybrid evolver interfaces and
      implementations, and hybrid orbit specialisations.
- [ ] Keep thin labelled adaptations of evolution-owned concepts in
      `evolution`, notably `LabelledEnclosure` and labelled orbit support.
- [ ] Keep `LabelledEnclosure` focused on labelling metadata and consistency
      around the numerical `Enclosure` operations; the numerical evolution
      logic remains in `Enclosure`.
- [ ] Use the `labelled` representation API exposed through `system -> hybrid`, without adding an unnecessary direct `evolution -> labelled` dependency.
- [ ] Replace System-owned `EvolverType`/`EnclosureType` associations with
      evolution-owned traits or equivalent dispatch, e.g. an
      `EvolutionTraits<SYS>` mapping.
- [ ] Verify that continuous, discrete, and hybrid evolvers consume System
      through its public system-description API only.

## H4 - Extract Analysis

`analysis` owns global reachability/safety algorithms and their global set
representations.

- [ ] Move `Storage`, reachability analysers, safety certificates, and their
      configurations into `analysis`.
- [ ] Move the corresponding hybrid analysis types into the same repository,
      including `HybridStorage` and `HybridReachabilityAnalyser`.
- [ ] Keep thin labelled adaptations of analysis-owned concepts in
      `analysis`, notably `LabelledStorage`.
- [ ] Keep `LabelledStorage` focused on state/auxiliary-space metadata,
      compatibility checks, labelled projections, and delegation to
      `Storage` for the actual paving operations.
- [ ] Replace System-owned `StorageType` associations with analysis-owned
      traits or equivalent dispatch, e.g. an `AnalysisTraits<SYS>` mapping.
- [ ] Keep state-space analysis aliases in analysis-owned traits rather than
      introducing them into the core `HybridSpace` representation.
- [ ] Keep `analysis -> evolution, smt, sat` as the high-level dependencies; use
      `labelled` transitively through `evolution -> system -> hybrid`.

## H5 - Classify source files by responsibility

Classify every file currently under `source/dynamics` and `source/hybrid`
by the concept it implements, not by whether its name contains `Hybrid` or
`Labelled`.

- [ ] Produce a file/class inventory for `hybrid`, `system`, `evolution`,
      and `analysis`.
- [ ] Retain the core hybrid representations, sets, paving, and graphics
      together in `source/hybrid`, analogously to `labelled`.
- [ ] Move hybrid evolution and analysis algorithms to their respective
      layers alongside continuous and discrete counterparts.
- [ ] Identify and split files that mix responsibilities before extracting
      any additional repositories.
- [ ] Remove obsolete forwarding headers only after consumers migrate.

The intended classification includes, at minimum:

```text
hybrid (representation layer):
    DiscreteLocation
    HybridSpace / MonolithicHybridSpace
    HybridConstraintSet / HybridBoundedConstraintSet
    Hybrid sets, paving, and graphics

system (system definitions):
    DiscreteEvent
    IteratedMap / VectorField / DifferentialInclusion
    HybridAutomaton / HybridAutomatonInterface

evolution (evolution of systems):
    Enclosure / LabelledEnclosure / HybridEnclosure
    EvolverInterface and concrete evolvers
    Simulators and Orbit specialisations

analysis (analysis of systems):
    Storage / LabelledStorage / HybridStorage
    ReachabilityAnalyser / HybridReachabilityAnalyser
    SafetyCertificate
```

## H6 - Enforce the high-level dependency contracts

- [ ] Add public-header isolation tests for `labelled`, `system`,
      `evolution`, and `analysis`.
- [ ] Add repository-level checks preventing `system` from depending on
      `evolution` or `analysis`; its lower-layer dependency on `hybrid`,
      including transitive `labelled`, is intentional.
- [ ] Verify that shared lower-level dependencies are reached transitively
      and consistently through the repository DAG.
- [ ] Verify that `analysis` receives `labelled` through `evolution` and
      does not add a direct dependency.
- [ ] Check all dependency diamonds against the same-commit consistency rule
      described in `structure.md`.
- [ ] Propagate each validated lower-level revision upward along every path
      before testing the aggregate `ariadne` repository.

## Cross-cutting packaging cleanup

- [ ] Remove Ariadne's remaining source-tree knowledge of Numeric used only to
      copy `FindGMP.cmake` and `FindMPFR.cmake`; keep the installed
      `AriadneConfig.cmake` contract valid while doing so.
- [ ] Continue using `ariadne_register_public_headers` and
      `ariadne_install_dependency_bundle` instead of manual transitive header
      installation lists.

## Validation discipline

For every boundary change:

1. compile affected public headers in isolation;
2. build the component without sibling include leakage;
3. consume only the immediate lower repository at the C++ boundary;
4. for Python, consume only the aggregate
   `pyariadne-<direct-dependency>` interface and utility header;
5. keep layer utility headers free of helpers already owned below;
6. run focused tests and standalone installation/consumer checks;
7. propagate the validated revision upward;
8. run Ariadne integration CI.

Do not treat successful aggregate Ariadne linkage as proof of a standalone
boundary.
