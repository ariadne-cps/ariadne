# Post-decoupling repository structure

This document describes the intended structure of Ariadne after the current decoupling effort. Repository names in prose and dependency descriptions use their actual lowercase, hyphenated spelling; the graph uses uppercase, space-separated display labels.

## Repository hierarchy

The top-level repository is `ariadne`. The list below and the diagram show a *transitively reduced* dependency graph: edges implied by other dependency paths are omitted. The graph therefore preserves transitive reachability without necessarily listing every direct dependency declaration in the repositories:

- `ariadne` uses `analysis`.
- `analysis` uses `evolution`, `smt`, `sat`. See [commit constraints](#analysis-commit-constraints).
- `evolution` uses `integration`, `system`, `labelling`, `symbolic-set`, `equation-solving`. See [commit constraints](#evolution-commit-constraints).
- `integration` uses `function`.
- `system` uses `symbolic`.
- `labelling` uses `symbolic`, `drawing`. See [commit constraints](#labelling-commit-constraints).
- `drawing` uses `rendering`, `paving`. See [commit constraints](#drawing-commit-constraints).
- `symbolic-set` uses `symbolic`, `set`. See [commit constraints](#symbolic-set-commit-constraints).
- `paving` uses `set`, `bdd`. See [commit constraints](#paving-commit-constraints).
- `set` uses `linear-programming`, `constraint-propagation`, `constraint-feasibility`. See [commit constraints](#set-commit-constraints).
- `smt` uses `symbolic`, `constraint-propagation`, `nonlinear-programming`, `equation-solving`, `threading`, `cdcl`. See [commit constraints](#smt-commit-constraints).
- `sat` uses `bdd`, `cdcl`. See [commit constraints](#sat-commit-constraints).
- `symbolic` uses `function`.
- `constraint-propagation` uses `function`.
- `constraint-feasibility` uses `nonlinear-programming`.
- `equation-solving` uses `function`.
- `nonlinear-programming` uses `function`.
- `linear-programming` uses `algebra`.
- `function` uses `geometry`.
- `geometry` uses `algebra`.
- `algebra` uses `interval`.
- `interval` uses `numeric`.
- `numeric` uses `foundation`, `logging`. See [commit constraints](#numeric-commit-constraints).
- `rendering` uses `utility`, `python-common`.
- `bdd` uses `threading`, `python-common`.
- `cdcl` uses `utility`, `python-common`.
- `threading` uses `logging`.
- `logging` uses `utility`.
- `foundation` uses `utility`, `python-common`.
- `utility` uses `macro`.

`python-common` and `macro` have no further dependencies in the graph. Each repository also uses `configuration`, which is omitted from the graph to avoid drawing an edge from every node.

All dependencies are exposed transitively. For example, `analysis` obtains `labelling` through `evolution`, so it does not need a direct `labelling` dependency. `symbolic-set` uses both `symbolic` and `set`; unlike the previous structure, `symbolic` uses `function` directly. `set` exposes `linear-programming`, `constraint-propagation`, and `constraint-feasibility` directly in the reduced graph, and `nonlinear-programming` transitively through `constraint-feasibility`. The repositories `constraint-propagation`, `equation-solving`, and `nonlinear-programming` use `function`, while `constraint-feasibility` uses `nonlinear-programming`.

## Commit consistency verification

The dependency graph is a directed acyclic graph (DAG). Each repository must use the same commit for any dependency exposed by different direct submodules, even if those submodules reach it through different paths. The checks below are derived from the **transitively reduced graph above** and assume that each submodule has already validated its own dependencies. The graph itself is unchanged.

### Verification rules

1. For each repository, consider every dependency reachable through two or more of its direct submodules. Only such shared dependencies can require a local cross-submodule commit comparison.
2. Compare commit identities (not merely repository names, branches, or version labels). A successful equality check on a dependency also fixes the commits of all its transitive dependencies, provided that the dependency's commit identifies an immutable, consistent submodule structure.
3. Perform the checks in dependency order: check an upstream/shared repository before using its equality to infer equality for repositories farther downstream. The local consistency of each submodule must also have been established before these inferences are used.
4. After checking a shared dependency, treat the direct submodules synchronized by that comparison as one equivalence group **for descendants of that shared dependency**. This does not assert that their unrelated dependencies match.
5. For a later shared dependency, compare just **one representative path from each still-distinct group** that reaches it. Choose the shortest path from any direct submodule in each group; if several paths are equally short, any one is sufficient.
6. Omit a downstream check if all direct submodules reaching that dependency already belong to one synchronized group. Do not omit a check that is still needed to connect another, previously unverified group.
7. `configuration`, which every repository uses but which is omitted from the graph, is excluded from the following graph-derived checks and must be handled separately.

All paths below start at a *direct submodule* of the repository named in the subsection. For each shared dependency, the listed paths must resolve to the same commit. The order of the subsections' checks is significant wherever a later check relies on an earlier equivalence.

**Example:** For `analysis`, checking `bdd` already synchronizes the `evolution` and `sat` branches for descendants of `bdd`. Therefore, the `threading` check only compares `sat -> bdd -> threading` with `smt -> threading`; a third `evolution` path is redundant.

### `analysis` commit constraints

- **`equation-solving`** (equal commits):
  - `evolution -> equation-solving`
  - `smt -> equation-solving`
- **`cdcl`** (equal commits):
  - `smt -> cdcl`
  - `sat -> cdcl`
- **`symbolic`** (equal commits):
  - `evolution -> system -> symbolic`
  - `smt -> symbolic`
- **`bdd`** (equal commits):
  - `evolution -> labelling -> drawing -> paving -> bdd`
  - `sat -> bdd`
- **`constraint-propagation`** (equal commits):
  - `evolution -> symbolic-set -> set -> constraint-propagation`
  - `smt -> constraint-propagation`
- **`threading`** (equal commits):
  - `sat -> bdd -> threading`
  - `smt -> threading`
- **`nonlinear-programming`** (equal commits):
  - `evolution -> symbolic-set -> set -> constraint-feasibility -> nonlinear-programming`
  - `smt -> nonlinear-programming`

### `evolution` commit constraints

- **`symbolic`** (equal commits):
  - `system -> symbolic`
  - `labelling -> symbolic`
  - `symbolic-set -> symbolic`
- **`set`** (equal commits):
  - `labelling -> drawing -> paving -> set`
  - `symbolic-set -> set`
- **`function`** (equal commits):
  - `integration -> function`
  - `system -> symbolic -> function`
  - `equation-solving -> function`

### `labelling` commit constraints

- **`function`** (equal commits):
  - `symbolic -> function`
  - `drawing -> paving -> set -> constraint-propagation -> function`

### `drawing` commit constraints

- **`python-common`** (equal commits):
  - `rendering -> python-common`
  - `paving -> bdd -> python-common`
- **`utility`** (equal commits):
  - `rendering -> utility`
  - `paving -> bdd -> threading -> logging -> utility`

### `symbolic-set` commit constraints

- **`function`** (equal commits):
  - `symbolic -> function`
  - `set -> constraint-propagation -> function`

### `paving` commit constraints

- **`logging`** (equal commits):
  - `set -> linear-programming -> algebra -> interval -> numeric -> logging`
  - `bdd -> threading -> logging`
- **`python-common`** (equal commits):
  - `set -> linear-programming -> algebra -> interval -> numeric -> foundation -> python-common`
  - `bdd -> python-common`

### `set` commit constraints

- **`function`** (equal commits):
  - `constraint-propagation -> function`
  - `constraint-feasibility -> nonlinear-programming -> function`
- **`algebra`** (equal commits):
  - `linear-programming -> algebra`
  - `constraint-propagation -> function -> geometry -> algebra`

### `smt` commit constraints

- **`function`** (equal commits):
  - `symbolic -> function`
  - `constraint-propagation -> function`
  - `nonlinear-programming -> function`
  - `equation-solving -> function`
- **`logging`** (equal commits):
  - `symbolic -> function -> geometry -> algebra -> interval -> numeric -> logging`
  - `threading -> logging`
- **`python-common`** (equal commits):
  - `symbolic -> function -> geometry -> algebra -> interval -> numeric -> foundation -> python-common`
  - `cdcl -> python-common`
- **`utility`** (equal commits):
  - `threading -> logging -> utility`
  - `cdcl -> utility`

### `sat` commit constraints

- **`python-common`** (equal commits):
  - `bdd -> python-common`
  - `cdcl -> python-common`
- **`utility`** (equal commits):
  - `bdd -> threading -> logging -> utility`
  - `cdcl -> utility`

### `numeric` commit constraints

- **`utility`** (equal commits):
  - `foundation -> utility`
  - `logging -> utility`

Repositories without cross-submodule constraints in the graph: `ariadne`, `integration`, `system`, `symbolic`, `constraint-propagation`, `constraint-feasibility`, `equation-solving`, `nonlinear-programming`, `linear-programming`, `function`, `geometry`, `algebra`, `interval`, `rendering`, `bdd`, `cdcl`, `threading`, `logging`, `foundation`, `utility`, `python-common`, `macro`.

In total, the graph yields **25 independent commit checks across 10 repositories**, using **54 representative paths**. The omitted downstream comparisons are implied by the verified commits and the dependency DAG.

## Dependency graph

```mermaid
graph TD
    ariadne["ARIADNE"]
    analysis["ANALYSIS"]
    evolution["EVOLUTION"]
    integration["INTEGRATION"]
    system["SYSTEM"]
    labelling["LABELLING"]
    drawing["DRAWING"]
    symbolic_set["SYMBOLIC SET"]
    paving["PAVING"]
    set["SET"]
    smt["SMT"]
    sat["SAT"]
    symbolic["SYMBOLIC"]
    constraint_propagation["CONSTRAINT PROPAGATION"]
    constraint_feasibility["CONSTRAINT FEASIBILITY"]
    equation_solving["EQUATION SOLVING"]
    nonlinear_programming["NONLINEAR PROGRAMMING"]
    linear_programming["LINEAR PROGRAMMING"]
    function["FUNCTION"]
    geometry["GEOMETRY"]
    algebra["ALGEBRA"]
    interval["INTERVAL"]
    numeric["NUMERIC"]
    rendering["RENDERING"]
    bdd["BDD"]
    cdcl["CDCL"]
    threading["THREADING"]
    logging["LOGGING"]
    foundation["FOUNDATION"]
    utility["UTILITY"]
    python_common["PYTHON COMMON"]
    macro["MACRO"]

    ariadne --> analysis
    analysis --> evolution
    analysis --> smt
    analysis --> sat
    evolution --> integration
    evolution --> system
    evolution --> labelling
    evolution --> symbolic_set
    evolution --> equation_solving
    integration --> function
    system --> symbolic
    labelling --> symbolic
    labelling --> drawing
    drawing --> rendering
    drawing --> paving
    symbolic_set --> symbolic
    symbolic_set --> set
    paving --> set
    paving --> bdd
    set --> linear_programming
    set --> constraint_propagation
    set --> constraint_feasibility
    smt --> symbolic
    smt --> constraint_propagation
    smt --> nonlinear_programming
    smt --> equation_solving
    smt --> threading
    smt --> cdcl
    sat --> bdd
    sat --> cdcl
    symbolic --> function
    constraint_propagation --> function
    constraint_feasibility --> nonlinear_programming
    equation_solving --> function
    nonlinear_programming --> function
    linear_programming --> algebra
    function --> geometry
    geometry --> algebra
    algebra --> interval
    interval --> numeric
    numeric --> foundation
    numeric --> logging
    rendering --> utility
    rendering --> python_common
    bdd --> threading
    bdd --> python_common
    cdcl --> utility
    cdcl --> python_common
    threading --> logging
    logging --> utility
    foundation --> utility
    foundation --> python_common
    utility --> macro
```

## Exposed dependencies and diamonds

Each repository exposes the repositories it uses to its consumers. Consequently, the dependency graph contains diamonds: the same lower-level repository may be reachable from a consumer through more than one dependency path.

These diamonds are intentional. At CMake level, consistency is enforced by checking that whenever the same dependency appears multiple times along dependency chains, every occurrence resolves to the same commit. A build must therefore not combine different commits of the same repository through different paths.

For example, `analysis` reaches `symbolic` through `evolution -> system`, `evolution -> labelling`, `evolution -> symbolic-set`, and `smt`. It reaches `set` through `evolution -> symbolic-set` and through `evolution -> labelling -> drawing -> paving`. The `set` dependencies introduce further diamonds for `constraint-propagation` and `constraint-feasibility`: `evolution` reaches `set` through both `symbolic-set` and `labelling -> drawing -> paving`. Additionally, `analysis` reaches `constraint-propagation` through `smt`.

The `nonlinear-programming` repository is reached through `set -> constraint-feasibility -> nonlinear-programming` and directly through `smt -> nonlinear-programming`. Similarly, `equation-solving` is reached directly from `evolution` and through `smt`. The shared `function` repository is exposed through `integration`, `symbolic`, `constraint-propagation`, `nonlinear-programming`, and `equation-solving`. Every such path must agree on the dependency commit.

The `sat` branch also creates a diamond for `bdd`: `analysis` reaches `bdd` through `sat` and through `evolution -> labelling -> drawing -> paving`. Moreover, `analysis` reaches `cdcl` through both `smt` and `sat`, and reaches `threading` directly through `smt` and transitively through `bdd`. All paths must agree on the commit of each shared repository.

`algebra` is reached through both `function -> geometry` and `set -> linear-programming`; `python-common` is shared by `rendering`, `bdd`, `cdcl`, and `foundation`.

## Minimal checkout strategy

The transitive dependency graph does not imply a recursive checkout of the whole submodule tree. The CMake dependency machinery checks out only what is strictly necessary for the repository being configured. When validating a dependency diamond, it descends into nested submodules only as far as needed to inspect and compare the referenced commits and verify that all paths resolve the same repository to the same commit.

Consequently, a recursive submodule checkout such as `git submodule update --init --recursive` is unnecessary for the normal build and consistency-checking workflow, and is inefficient because it materializes repositories that the current configuration does not need.

## Propagating changes upward

When a repository is modified and a repository above it must see those changes, the new dependency commit must be propagated upward along **every** path that reaches the consumer.

This is necessary because a repository can be exposed through several branches of the graph. Updating only one branch would make different paths resolve the same dependency to different commits, which the CMake consistency check is expected to reject.

Representative (not exhaustive) paths for a change to `function` to reach `ariadne` include:

```text
function -> integration                                             -> evolution -> analysis -> ariadne
function -> symbolic -> system                                     -> evolution -> analysis -> ariadne
function -> symbolic -> labelling                                  -> evolution -> analysis -> ariadne
function -> symbolic -> symbolic-set                               -> evolution -> analysis -> ariadne
function -> constraint-propagation -> set -> symbolic-set           -> evolution -> analysis -> ariadne
function -> nonlinear-programming -> constraint-feasibility -> set -> symbolic-set -> evolution -> analysis -> ariadne
function -> equation-solving                                      -> evolution -> analysis -> ariadne
function -> symbolic -> smt                                                     -> analysis -> ariadne
function -> nonlinear-programming -> smt                                        -> analysis -> ariadne
```

Additional paths pass through `paving`, `set`, and their consumers. The representative routes above do not replace the requirement to propagate along **all** paths in the dependency graph.

A change to `bdd` must be propagated through both `sat` and `paving -> drawing -> labelling -> evolution`, including all upstream consumers. A change to `threading` must additionally propagate through `smt`, alongside the routes through `bdd`. A change to `cdcl` must propagate through both `smt` and `sat` before reaching `analysis`. Similarly, changes to `constraint-propagation`, `constraint-feasibility`, `equation-solving`, or `nonlinear-programming` require updating every direct and transitive consumer shown in the graph.

The special case is `configuration`: because every repository depends on it, making a new `configuration` commit visible at the top requires propagating that commit through the entire dependency graph, all the way to `ariadne`.
