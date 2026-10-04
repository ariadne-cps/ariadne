# Decoupling Ariadne components

This directory records the architectural state and the remaining work for
splitting Ariadne into independently buildable repositories.

Dependency arrows use **provider -> dependent**.

## Current state

The low-level repository chain through Algebra is complete and validated:

```text
utility
  -> foundation
       -> numeric
            -> interval
                 -> algebra
                      -> function   (next extraction)
```

`foundation` is the renamed former `paradigm` repository. The semantic C++
concepts `Paradigm`, `ParadigmCode`, `ParadigmTraits`, etc. retain their
names.

Ariadne consumes standalone Algebra rather than owning a local Algebra
implementation. Function is the next component to be decoupled and extracted.

See:

- [Current state](current-state.md) for the active architecture and boundary
  rules.
- [Working plan](work-plan.md) for the remaining Function extraction work.

## Historical baseline

The following files describe the repository baseline captured on 2026-09-29 at
`949d7e044ae65837fc02e6387701b10e7c15ddc6`. They are retained as historical
evidence and are not the current dependency graph:

- [Coupling analysis](coupling-analysis.md)
- [Dependency graph](dependency-graph.svg)
- [Dependency summary](dependency-summary.csv)
- [Include evidence](include-evidence.csv)

The baseline data is intentionally not rewritten as the architecture changes.

## Boundary rules

Each repository owns only its layer and consumes its immediate lower repository.
The same rule applies to Python bindings: a layer consumes the aggregate
`pyariadne-<component>` interface and utility header of its direct dependency,
without reaching into nested submodules.

Repositories that use shared CMake infrastructure carry
`submodules/configuration` directly. When the same Configuration repository is
also present transitively, `require_same_dependency_commit` verifies that the
gitlinks agree.

Public headers are registered with `ariadne_register_public_headers` and
installed transitively with `ariadne_install_dependency_bundle`. Consumers
should not reproduce nested header-install lists manually.
