# Post-decoupling repository structure

This document describes the intended structure of Ariadne after the current decoupling effort.

## Repository hierarchy

The top-level repository is `ariadne`. Its dependency chain is organized as follows:

- `ariadne` uses `analysis`.
- `analysis` uses both `evolution` and `satisfiability`.
- `evolution` uses `integration`, `system`, and `labelling`.
- `integration` uses `kernel`.
- `system` uses `symbolic`, but does not use `labelling` directly.
- `satisfiability` uses `symbolic`.
- `labelling` uses `symbolic`.
- `symbolic` uses `kernel`.
- `kernel` remains an aggregate repository internally composed of the modules `function`, `geometry`, `drawing`, and `solving`. Until these modules can be decoupled further, they remain part of `kernel`.
- `kernel` uses `bdd` and `algebra`.
- `bdd` uses `threading` and `python-common`.
- `threading` uses `logging`.
- `logging` uses `utility`.
- `utility` uses `macro`.
- `algebra` uses `interval`.
- `interval` uses `numeric`.
- `numeric` uses `foundation` and `logging`.
- `foundation` uses `utility` and `python-common`.

In addition, every repository uses `configuration`. This dependency is omitted from the graph below because it would otherwise appear from every repository node.

Dependencies are exposed transitively. In particular, `analysis` does not need a direct dependency on `labelling`: it receives `labelling` through `evolution`. Likewise, `system`, `labelling`, and `satisfiability` do not need direct dependencies on `kernel`, since they receive it through `symbolic`.

## Dependency graph

In the graph, a solid arrow means "uses", while a dotted arrow marks an internal module contained in the `kernel` aggregate.

```mermaid
graph TD
    ariadne --> analysis

    analysis --> evolution
    analysis --> satisfiability

    evolution --> integration
    evolution --> system
    evolution --> labelling

    integration --> kernel

    system --> symbolic

    satisfiability --> symbolic

    labelling --> symbolic

    symbolic --> kernel

    kernel --> bdd
    kernel --> algebra

    kernel -. internal module .-> function
    kernel -. internal module .-> geometry
    kernel -. internal module .-> drawing
    kernel -. internal module .-> solving

    bdd --> threading
    bdd --> python_common["python-common"]

    threading --> logging
    logging --> utility
    utility --> macro

    algebra --> interval
    interval --> numeric
    numeric --> foundation
    numeric --> logging

    foundation --> utility
    foundation --> python_common
```

## Exposed dependencies and diamonds

Each repository exposes the repositories it uses to its consumers. Consequently, the dependency graph can contain diamonds: the same lower-level repository may be reachable from a consumer through more than one dependency path.

These diamonds are intentional. At CMake level, consistency is enforced by checking that whenever the same dependency appears multiple times along the dependency chains, every occurrence resolves to the same commit. A build must therefore not combine different commits of the same repository through different paths.

For example, `kernel` is reached from `evolution` directly through `integration`, and transitively through `symbolic` along both the `system` and `labelling` branches. It is also reached from `analysis` through `satisfiability -> symbolic`. All paths must ultimately expose the same `kernel` commit.

Similarly, `analysis` reaches `labelling` transitively through `evolution`; no additional direct `analysis -> labelling` dependency is required.

## Minimal checkout strategy

The transitive dependency graph does not imply a recursive checkout of the whole submodule tree. The CMake dependency machinery checks out only what is strictly necessary for the repository being configured. When validating a dependency diamond, it descends into nested submodules only as far as needed to inspect and compare the referenced commits and verify that all paths resolve the same repository to the same commit.

Consequently, a recursive submodule checkout such as `git submodule update --init --recursive` is unnecessary for the normal build and consistency-checking workflow, and is inefficient because it materializes repositories that the current configuration does not need.

## Propagating changes upward

When a repository is modified and a repository above it must see those changes, the new dependency commit must be propagated upward along **every** path that reaches the consumer.

This is necessary because a repository can be exposed through several branches of the graph. Updating only one branch would make different paths resolve the same dependency to different commits, which the CMake consistency check is expected to reject.

As a concrete example, a change to `kernel` that must become visible to `ariadne` has to be propagated through all paths that reach it:

```text
kernel -> integration                 -> evolution -> analysis -> ariadne
kernel -> symbolic -> system          -> evolution -> analysis -> ariadne
kernel -> symbolic -> labelling        -> evolution -> analysis -> ariadne
kernel -> symbolic -> satisfiability              -> analysis -> ariadne
```

The special case is `configuration`: because every repository depends on it, making a new `configuration` commit visible at the top requires propagating that commit through the entire dependency graph, all the way to `ariadne`.
