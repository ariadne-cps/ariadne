# Structure

Repository names in prose and dependency descriptions use their actual lowercase, hyphenated spelling; the graph uses uppercase, space-separated display labels.

```mermaid
%%{init: {"flowchart": {"nodeSpacing": 8, "rankSpacing": 90}}}%%
flowchart TD
    ariadne["ARIADNE"]
    reasoning["REASONING"]
    sat["SAT"]
    smt["SMT"]
    neural_set["NEURAL<br/>SET"]
    evolution["EVOLUTION"]
    integration["INTEGRATION"]
    system["SYSTEM"]
    hybrid["HYBRID"]
    labelled["LABELLED"]
    drawing["DRAWING"]
    expression_set["EXPRESSION<br/>SET"]
    paving["PAVING"]
    tree["TREE"]
    set["SET"]
    ann["ANN"]
    symbolic["SYMBOLIC"]
    constraint_solving["CONSTRAINT<br/>SOLVING"]
    equation_solving["EQUATION<br/>SOLVING"]
    nonlinear_programming["NONLINEAR<br/>PROGRAMMING"]
    linear_programming["LINEAR<br/>PROGRAMMING"]
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
    python_common["PYTHON<br/>COMMON"]

    ariadne --> reasoning
    reasoning --> sat
    reasoning --> smt
    reasoning --> evolution

    %% Layout-only links: these invisible edges adjust vertical ranks,
    %% but do not represent repository dependencies.
    sat ~~~ smt
    smt ~~~ evolution
    equation_solving ~~~ integration
    symbolic ~~~ paving
    tree ~~~ bdd
    constraint_solving ~~~ threading
    sat --> bdd
    sat --> cdcl
    smt --> symbolic
    smt --> cdcl
    smt --> constraint_solving
    smt --> equation_solving

    neural_set --> ann
    neural_set --> set
    ann --> threading

    evolution --> integration
    evolution --> system
    evolution --> equation_solving
    integration --> function
    system --> hybrid
    system --> neural_set
    hybrid --> labelled
    labelled --> expression_set
    labelled --> drawing

    expression_set --> symbolic
    expression_set --> set
    drawing --> rendering
    drawing --> paving
    paving --> set
    paving --> bdd
    paving --> tree
    tree --> foundation
    tree --> threading
    set --> linear_programming
    set --> constraint_solving

    symbolic --> function
    constraint_solving --> nonlinear_programming
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
    cdcl --> threading
    cdcl --> python_common
    threading --> logging
    logging --> utility
    foundation --> utility
    foundation --> python_common
```

The list below and the diagram show a *transitively reduced* dependency graph: edges implied by other dependency paths are omitted. The graph therefore preserves transitive reachability without necessarily listing every direct dependency declaration in the repositories:

- `ariadne` uses `reasoning`.
- `reasoning` uses `evolution`, `smt`, `sat`. See [commit constraints](#reasoning-commit-constraints).
- `neural-set` uses `ann`, `set`. See [commit constraints](#neural-set-commit-constraints).
- `ann` uses `threading`.
- `evolution` uses `integration`, `system`, `equation-solving`. See [commit constraints](#evolution-commit-constraints).
- `integration` uses `function`.
- `system` uses `hybrid`, `neural-set`. See [commit constraints](#system-commit-constraints).
- `hybrid` uses `labelled`.
- `labelled` uses `expression-set`, `drawing`. See [commit constraints](#labelled-commit-constraints).
- `drawing` uses `rendering`, `paving`. See [commit constraints](#drawing-commit-constraints).
- `expression-set` uses `symbolic`, `set`. See [commit constraints](#expression-set-commit-constraints).
- `paving` uses `set`, `bdd`, `tree`. See [commit constraints](#paving-commit-constraints).
- `tree` uses `foundation`, `threading`.
- `set` uses `linear-programming`, `constraint-solving`. See [commit constraints](#set-commit-constraints).
- `smt` uses `symbolic`, `constraint-solving`, `equation-solving`, `cdcl`. See [commit constraints](#smt-commit-constraints).
- `sat` uses `bdd`, `cdcl`. See [commit constraints](#sat-commit-constraints).
- `symbolic` uses `function`.
- `constraint-solving` uses `nonlinear-programming`.
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
- `cdcl` uses `threading`, `python-common`.
- `threading` uses `logging`.
- `logging` uses `utility`.
- `foundation` uses `utility`, `python-common`.

`python-common` has no further dependencies in the graph. Each repository also uses `build`, which is omitted from the graph to avoid drawing an edge from every node.

The `ariadne-cps/testing` repository is included as a test-only Git submodule in repositories that define C++ test suites. The top-level `ariadne` aggregator has no separate test executables and obtains Testing transitively through `reasoning`. Testing is therefore omitted from the functional dependency graph and graph-derived commit consistency checks.

All dependencies are exposed transitively. For example, `reasoning` obtains `labelled` through `evolution -> system -> hybrid`, so it does not need a direct `labelled` dependency. Similarly, `evolution` reaches `labelled` through `system -> hybrid`, and `system` reaches `symbolic` through `hybrid -> labelled -> expression-set`; none of these transitive relationships needs an additional edge in the reduced graph. These reduced edges may still be declared directly by individual repositories; the graph represents reachability, not every direct declaration.

## Commit consistency verification

The dependency graph is a directed acyclic graph (DAG). Each repository must use the same commit for any dependency exposed by different direct submodules, even if those submodules reach it through different paths. The checks below are derived from the **transitively reduced graph above** and assume that each submodule has already validated its own dependencies. The checks below reflect the graph shown above.

### Verification rules

1. For each repository, consider every dependency reachable through two or more of its direct submodules. Only such shared dependencies can require a local cross-submodule commit comparison.
2. Compare commit identities (not merely repository names, branches, or version labels). A successful equality check on a dependency also fixes the commits of all its transitive dependencies, provided that the dependency's commit identifies an immutable, consistent submodule structure.
3. Perform the checks in dependency order: check an upstream/shared repository before using its equality to infer equality for repositories farther downstream. The local consistency of each submodule must also have been established before these inferences are used.
4. After checking a shared dependency, treat the direct submodules synchronized by that comparison as one equivalence group **for descendants of that shared dependency**. This does not assert that their unrelated dependencies match.
5. For a later shared dependency, compare just **one representative path from each still-distinct group** that reaches it. Choose the shortest path from any direct submodule in each group; if several paths are equally short, any one is sufficient.
6. Omit a downstream check if all direct submodules reaching that dependency already belong to one synchronized group. Do not omit a check that is still needed to connect another, previously unverified group.
7. `build`, which every repository uses but which is omitted from the graph, is excluded from the following graph-derived checks and must be handled separately.

All paths below start at a *direct submodule* of the repository named in the subsection. For each shared dependency, the listed paths must resolve to the same commit. The order of the subsections' checks is significant wherever a later check relies on an earlier equivalence.

**Example:** For `reasoning`, the `bdd` check synchronizes the `evolution` and `sat` paths to `bdd`, and the `cdcl` check synchronizes the `smt` and `sat` paths to `cdcl`. Since `sat` has already verified that `bdd -> threading` and `cdcl -> threading` refer to the same commit, these three submodules already agree on `threading`, so `reasoning` needs no separate `threading` check.

### `reasoning` commit constraints

- **`equation-solving`**:
  - `evolution -> equation-solving`
  - `smt -> equation-solving`
- **`cdcl`**:
  - `smt -> cdcl`
  - `sat -> cdcl`
- **`symbolic`**:
  - `evolution -> system -> hybrid -> labelled -> expression-set -> symbolic`
  - `smt -> symbolic`
- **`bdd`**:
  - `evolution -> system -> hybrid -> labelled -> drawing -> paving -> bdd`
  - `sat -> bdd`
- **`constraint-solving`**:
  - `evolution -> system -> hybrid -> labelled -> expression-set -> set -> constraint-solving`
  - `smt -> constraint-solving`

### `neural-set` commit constraints

- **`logging`**:
  - `ann -> threading -> logging`
  - `set -> linear-programming -> algebra -> interval -> numeric -> logging`

### `evolution` commit constraints

- **`function`**:
  - `integration -> function`
  - `system -> hybrid -> labelled -> expression-set -> symbolic -> function`
  - `equation-solving -> function`

### `system` commit constraints

- **`set`**:
  - `hybrid -> labelled -> expression-set -> set`
  - `neural-set -> set`
- **`threading`**:
  - `hybrid -> labelled -> drawing -> paving -> bdd -> threading`
  - `neural-set -> ann -> threading`

### `labelled` commit constraints

- **`set`**:
  - `expression-set -> set`
  - `drawing -> paving -> set`

### `drawing` commit constraints

- **`python-common`**:
  - `rendering -> python-common`
  - `paving -> bdd -> python-common`
- **`utility`**:
  - `rendering -> utility`
  - `paving -> tree -> foundation -> utility`

### `expression-set` commit constraints

- **`function`**:
  - `symbolic -> function`
  - `set -> constraint-solving -> nonlinear-programming -> function`

### `paving` commit constraints

- **`foundation`**:
  - `set -> linear-programming -> algebra -> interval -> numeric -> foundation`
  - `tree -> foundation`
- **`threading`**:
  - `bdd -> threading`
  - `tree -> threading`
- **`logging`**:
  - `set -> linear-programming -> algebra -> interval -> numeric -> logging`
  - `bdd -> threading -> logging`
- **`python-common`**:
  - `tree -> foundation -> python-common`
  - `bdd -> python-common`

### `set` commit constraints

- **`algebra`**:
  - `linear-programming -> algebra`
  - `constraint-solving -> nonlinear-programming -> function -> geometry -> algebra`

### `smt` commit constraints

- **`function`**:
  - `symbolic -> function`
  - `constraint-solving -> nonlinear-programming -> function`
  - `equation-solving -> function`
- **`logging`**:
  - `symbolic -> function -> geometry -> algebra -> interval -> numeric -> logging`
  - `cdcl -> threading -> logging`
- **`python-common`**:
  - `symbolic -> function -> geometry -> algebra -> interval -> numeric -> foundation -> python-common`
  - `cdcl -> python-common`
### `sat` commit constraints

- **`threading`**:
  - `bdd -> threading`
  - `cdcl -> threading`
- **`python-common`**:
  - `bdd -> python-common`
  - `cdcl -> python-common`

### `numeric` commit constraints

- **`utility`**:
  - `foundation -> utility`
  - `logging -> utility`

## Minimal checkout strategy

The transitive dependency graph does not imply a recursive checkout of the whole submodule tree. The CMake dependency machinery checks out only what is strictly necessary for the repository being configured. When validating a dependency diamond, it descends into nested submodules only as far as needed to inspect and compare the referenced commits and verify that all paths resolve the same repository to the same commit.

Consequently, a recursive submodule checkout such as `git submodule update --init --recursive` is unnecessary for the normal build and consistency-checking workflow, and is inefficient because it materializes repositories that the current configuration does not need.

## Propagating changes upward

When a repository is modified and a repository above it must see those changes, the new dependency commit must be propagated upward along **every** path that reaches the consumer.

This is necessary because a repository can be exposed through several branches of the graph. Updating only one branch would make different paths resolve the same dependency to different commits, which the CMake consistency check is expected to reject.

The special case is `build`: because every repository depends on it, making a new `build` commit visible at the top requires propagating that commit through the entire dependency graph, all the way to `ariadne`.
