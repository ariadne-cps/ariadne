# Ariadne decoupling work plan

Status: active  
Base state: low-level chain through Algebra extracted and green  
Last updated: 2026-10-13

## Established chain

```text
utility -> foundation -> numeric -> interval -> algebra -> function
```

Utility, Foundation, Numeric, Interval and Algebra are standalone. The completed
Foundation/Algebra extraction history is recorded in Git history and
`current-state.md`; completed migration checklists are intentionally omitted
from this active plan.

## Target repository structure

The current target dependency graph is documented in `structure.md`. At the
high level, the intended direct dependencies are:

```text
symbolic       -> kernel
labelled      -> symbolic
system         -> symbolic
integration    -> kernel
evolution      -> integration, system, labelled
satisfiability -> symbolic
analysis       -> evolution, satisfiability
ariadne        -> analysis
```

Dependencies are exposed transitively. In particular:

- `system`, `labelled`, and `satisfiability` obtain `kernel` through
  `symbolic`; they do not declare a direct `kernel` dependency.
- `analysis` obtains `labelled` through `evolution`; it does not declare a
  direct `labelled` dependency.
- `system` must not depend on `labelled`, `evolution`, or `analysis`.

The high-level split is by responsibility rather than by the current
`dynamics/` and `hybrid/` source directories:

- `system` describes discrete maps, vector fields, differential inclusions,
  hybrid automata, locations, events, and state spaces;
- `evolution` contains enclosures, evolvers, simulators, and orbit types;
- `analysis` contains global storage, reachability analysis, safety
  certificates, and related algorithms.

Hybrid types follow the same responsibility split instead of forming a
standalone repository.

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

## H2 - Extract System without Labelled or upward dependencies

The `system` repository must describe what evolves without knowing how it is
evolved or analysed.

- [ ] Move the system-description classes from the current Dynamics and Hybrid
      modules into `system`: `IteratedMap`, `VectorField`,
      `DifferentialInclusion`, `HybridAutomaton`,
      `HybridAutomatonInterface`, `DiscreteLocation`, `DiscreteEvent`,
      `HybridSpace`, and directly related specification types.
- [ ] Remove `EvolverType`, `EnclosureType`, and `StorageType` aliases from
      `IteratedMap`, `VectorField`, `DifferentialInclusion`, and
      `HybridAutomatonInterface`.
- [ ] Remove corresponding forward declarations and includes of evolver,
      enclosure, storage, and `Labelled*` types from System public headers.
- [ ] Remove the unused `LabelledSet` forward declaration from
      `vector_field.hpp`.
- [ ] Keep system state-variable semantics expressed through `symbolic`
      concepts such as `RealSpace`, variables, assignments, and predicates.
- [ ] Keep `DifferentialInclusion` input domains expressible with
      `RealVariablesBox` without introducing a `labelled` dependency.
- [ ] Reduce `HybridSpace` to the system-level notion of locations mapped to
      `RealSpace`.
- [ ] Move analysis-oriented `HybridSpace` aliases such as
      `BoundingDomainType`, set-interface aliases, and
      `SetApproximationType` out of `system` and into analysis-side traits.
- [ ] Remove indirect includes such as `hybrid_set_interface.hpp` from
      `HybridSpace` once those analysis responsibilities have moved.
- [ ] Add a boundary check proving that System public headers do not include or
      require `labelled`, `evolution`, or `analysis`.

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
- [ ] Make the direct dependency `evolution -> labelled` explicit.
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
- [ ] Move state-space analysis aliases currently stored in `HybridSpace` to
      analysis-owned traits.
- [ ] Make `analysis -> evolution` and `analysis -> satisfiability` the
      direct high-level dependencies; use the exposed `labelled` dependency
      inherited through `evolution` instead of adding a direct repository
      edge.

## H5 - Reclassify current Dynamics and Hybrid source files

Before physically extracting repositories, classify every file currently under
`source/dynamics` and `source/hybrid` by its target responsibility.

- [ ] Produce a file/class inventory for `system`, `evolution`, and
      `analysis`.
- [ ] Treat "hybrid" as a variant within each responsibility rather than as a
      target repository boundary.
- [ ] Treat the `Labelled*` prefix as an implementation/API variant, not as a
      repository selector: the owner is determined by the underlying concept.
- [ ] Identify files that mix responsibilities and split them before moving
      repositories.
- [ ] Remove obsolete cross-module forwarding headers only after consumers have
      migrated.

The intended classification includes, at minimum:

```text
system:
    IteratedMap
    VectorField
    DifferentialInclusion
    HybridAutomaton / HybridAutomatonInterface
    DiscreteLocation / DiscreteEvent
    HybridSpace

evolution:
    Enclosure / LabelledEnclosure / HybridEnclosure
    EvolverInterface and concrete evolvers
    simulators
    Orbit specialisations

analysis:
    Storage / LabelledStorage / HybridStorage
    ReachabilityAnalyser / HybridReachabilityAnalyser
    SafetyCertificate
```

## H6 - Enforce the high-level dependency contracts

- [ ] Add public-header isolation tests for `labelled`, `system`,
      `evolution`, and `analysis`.
- [ ] Add repository-level checks preventing `system` from acquiring a
      `labelled`, `evolution`, or `analysis` dependency.
- [ ] Verify that `system`, `labelled`, and `satisfiability` see
      `kernel` only through the exposed `symbolic` dependency.
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
