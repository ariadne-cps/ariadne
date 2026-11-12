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
    ann --> python_common

    evolution --> integration
    evolution --> system
    evolution --> equation_solving
    integration --> function
    integration --> logging
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
    equation_solving --> logging
    nonlinear_programming --> function
    nonlinear_programming --> logging
    linear_programming --> algebra
    linear_programming --> logging
    function --> geometry
    geometry --> algebra
    algebra --> interval
    interval --> numeric
    numeric --> foundation

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

The list below and the diagram show a *transitively reduced* functional dependency graph: edges implied by other dependency paths are omitted. The graph therefore preserves transitive reachability without necessarily listing every direct dependency declaration in the repositories:

- `ariadne` uses `reasoning`.
- `reasoning` uses `evolution`, `smt`, `sat`. See [commit constraints](#reasoning-commit-constraints).
- `neural-set` uses `ann`, `set`. See [commit constraints](#neural-set-commit-constraints).
- `ann` uses `threading`, `python-common`.
- `evolution` uses `integration`, `system`, `equation-solving`. See [commit constraints](#evolution-commit-constraints).
- `integration` uses `function`, `logging`.
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
- `equation-solving` uses `function`, `logging`.
- `nonlinear-programming` uses `function`, `logging`.
- `linear-programming` uses `algebra`, `logging`.
- `function` uses `geometry`.
- `geometry` uses `algebra`.
- `algebra` uses `interval`.
- `interval` uses `numeric`.
- `numeric` uses `foundation`.
- `rendering` uses `utility`, `python-common`.
- `bdd` uses `threading`, `python-common`.
- `cdcl` uses `threading`, `python-common`.
- `threading` uses `logging`.
- `logging` uses `utility`.
- `foundation` uses `utility`, `python-common`.

`python-common` has no further dependencies in the graph. Every repository also uses `build`, which is omitted from the graph to avoid drawing an edge from every node.

`testing` depends only on `build`. Every repository except `python-common` and the top-level `ariadne` includes `testing`, but a compilation uses only the `testing` submodule of its root repository; copies nested in dependencies are not used. Consequently, the commits of different `testing` submodules do not need to match. The `build` commit referenced by the active `testing` submodule **must**, however, be consistent with all other `build` commits used in the compilation.

A selected subset of repositories also includes `benchmarking`. Its functional dependencies are `logging` and `build`; like other repositories except the two noted above, it also includes `testing`. Both `testing` and `benchmarking` are omitted from the functional dependency graph. Their implications for commit consistency are described below.

All dependencies are exposed transitively. For example, `reasoning` obtains `labelled` through `evolution -> system -> hybrid`, so it does not need a direct `labelled` dependency. Similarly, `evolution` reaches `labelled` through `system -> hybrid`, and `system` reaches `symbolic` through `hybrid -> labelled -> expression-set`; none of these transitive relationships needs an additional edge in the reduced graph. These reduced edges may still be declared directly by individual repositories; the graph represents reachability, not every direct declaration.

## Commit consistency verification

The dependency graph is a directed acyclic graph (DAG). Each repository must use the same commit for any dependency exposed by different direct submodules, even if those submodules reach it through different paths. The checks below are derived from the **transitively reduced functional graph above** and assume that each submodule has already validated its own dependencies. Auxiliary repositories omitted from the graph are handled by separate rules.

### Verification rules

1. For each repository, consider every dependency reachable through two or more of its direct submodules. Only such shared dependencies can require a local cross-submodule commit comparison.
2. Compare commit identities (not merely repository names, branches, or version labels). A successful equality check on a dependency also fixes the commits of all its transitive dependencies, provided that the dependency's commit identifies an immutable, consistent submodule structure.
3. Perform the checks in dependency order: check an upstream/shared repository before using its equality to infer equality for repositories farther downstream. The local consistency of each submodule must also have been established before these inferences are used.
4. After checking a shared dependency, treat the direct submodules synchronized by that comparison as one equivalence group **for descendants of that shared dependency**. This does not assert that their unrelated dependencies match.
5. For a later shared dependency, compare just **one representative path from each still-distinct group** that reaches it. Choose the shortest path from any direct submodule in each group; if several paths are equally short, any one is sufficient.
6. Omit a downstream check if all direct submodules reaching that dependency already belong to one synchronized group. Do not omit a check that is still needed to connect another, previously unverified group.
7. `build`, which every repository uses but which is omitted from the graph, must be checked separately. In particular, the comparison must include `testing -> build` for the active `testing` submodule of the compilation root and, wherever `benchmarking` is included, `benchmarking -> build`, alongside all other paths to `build`.
8. No consistency check between distinct **`testing` commits** is required: only the compilation root's `testing` is used. This exemption does not apply to its `build` dependency, which must pass rule 7.
9. Wherever `benchmarking` is included, its `logging` dependency also requires consistency: add `benchmarking -> logging` to the `logging` constraint if the root repository uses `logging`; otherwise add `benchmarking -> logging -> utility` to the `utility` constraint, since the root necessarily uses `utility`. These comparisons are in addition to the functional graph-derived constraints.

The following per-repository paths cover the functional graph; the auxiliary checks in rules 7 and 9 must also be applied where relevant. All paths below start at a *direct submodule* of the repository named in the subsection. For each shared dependency, the listed paths must resolve to the same commit. The order of the subsections' checks is significant wherever a later check relies on an earlier equivalence.

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
  - `set -> linear-programming -> logging`
- **`python-common`**:
  - `ann -> python-common`
  - `set -> linear-programming -> algebra -> interval -> numeric -> foundation -> python-common`

### `evolution` commit constraints

- **`function`**:
  - `integration -> function`
  - `system -> hybrid -> labelled -> expression-set -> symbolic -> function`
  - `equation-solving -> function`
- **`logging`**:
  - `integration -> logging`
  - `system -> neural-set -> ann -> threading -> logging`
  - `equation-solving -> logging`

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
  - `set -> linear-programming -> logging`
  - `bdd -> threading -> logging`
- **`python-common`**:
  - `tree -> foundation -> python-common`
  - `bdd -> python-common`

### `set` commit constraints

- **`algebra`**:
  - `linear-programming -> algebra`
  - `constraint-solving -> nonlinear-programming -> function -> geometry -> algebra`
- **`logging`**:
  - `linear-programming -> logging`
  - `constraint-solving -> nonlinear-programming -> logging`

### `smt` commit constraints

- **`function`**:
  - `symbolic -> function`
  - `constraint-solving -> nonlinear-programming -> function`
  - `equation-solving -> function`
- **`logging`**:
  - `constraint-solving -> nonlinear-programming -> logging`
  - `equation-solving -> logging`
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

## Minimal checkout strategy

The transitive dependency graph does not imply a recursive checkout of the whole submodule tree. The CMake dependency machinery checks out only what is strictly necessary for the repository being configured. When validating a dependency diamond, it descends into nested submodules only as far as needed to inspect and compare the referenced commits and verify that all paths resolve the same repository to the same commit.

**Each logical dependency is certified only once per CMake configuration run.** On the first successful `require_same_dependency_commit(NAME CHAINS ...)` call, CMake compares all the specified Gitlink commits and records a `GLOBAL` certification for `NAME`. Subsequent checks for the same logical repository in nested projects return immediately: they neither repeat the comparisons nor initialize submodules solely to follow those chains. Certification is recorded *only after* all paths agree; a mismatch still fails configuration. The property is keyed by the logical dependency name, **not by an individual Git commit SHA or checkout path**, and is reset for the next configuration run.

This relies on the repository invariant that each pinned submodule commit is already internally consistent. For example, once `paving` has certified `logging` across the paths `set -> linear-programming -> logging` and `bdd -> threading -> logging`, later nested checks for that same logical dependency can reuse the certification. This does not exempt the root's active `testing -> build` path or the relevant `benchmarking` paths from their required consistency checks, nor does it allow skipping an uncertified cross-branch comparison. Skipping a repeated consistency check does **not** suppress checkout of a submodule that must actually be configured or built.

Consequently, a recursive submodule checkout such as `git submodule update --init --recursive` is unnecessary for the normal build and consistency-checking workflow, and is inefficient because it materializes repositories that the current configuration does not need.

## Propagating changes upward

When a repository is modified and a repository above it must see those changes, the new dependency commit must be propagated upward along **every** path that reaches the consumer.

This is necessary because a repository can be exposed through several branches of the graph. Updating only one branch would make different paths resolve the same dependency to different commits, which the CMake consistency check is expected to reject.

The special case is `build`: because every repository depends on it, making a new `build` commit visible at the top requires propagating that commit through the entire dependency graph, all the way to `ariadne`. This also includes maintaining consistency with the `build` commits pinned by the compilation root's `testing` and by any included `benchmarking`; `benchmarking` must additionally remain consistent on `logging` or, where appropriate, `utility`.

### Coordinated work branches and fast-forward integration

When changes span `ariadne` and multiple dependency repositories, development may take place on corresponding work branches in each repository.

During development, newly required dependency commits must be propagated **bottom-up through the work branches**, updating all affected gitlinks along every relevant dependency path and validating commit consistency before integration.

If each work branch is subsequently integrated into its repository's `main` through a **fast-forward merge** (`git merge --ff-only`), the resulting `main` tip is exactly the former tip of the work branch. Since gitlinks record immutable commit SHAs rather than branch names, all dependency references remain unchanged.

**Consequently, no additional bottom-up gitlink propagation is required after fast-forward integration**, provided that the work branches already contain the intended, mutually consistent dependency commits.

A fast-forward does not automatically advance gitlinks to the latest `main` tips of their dependencies. If integration rewrites commits (e.g., through rebasing or squashing), introduces new dependency commits that must be referenced, or otherwise changes the intended dependency pins, the affected gitlinks must be updated and consistency revalidated. A non-fast-forward merge does not necessarily invalidate existing gitlinks, but the new `main` tip will no longer equal the merged branch tip.

When publishing the integration, push dependency repositories before their consumers, following the dependency graph bottom-up. This ensures that every referenced commit is already available remotely when the corresponding parent repository is published.
