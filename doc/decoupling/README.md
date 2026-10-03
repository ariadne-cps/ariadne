# Decoupling Ariadne components

This directory records the architectural analysis and current work for moving
low-level Ariadne components into independently buildable repositories.

Dependency arrows in the historical analysis use **provider -> dependent**: an
arrow points to the component that consumes the provider.

## Read this first

- [Current state](current-state.md): verified repository structure as of
  2026-10-03, the target Foundation boundary, and the residual
  `algebra`/`function` dependencies that must be removed.
- [Working plan](work-plan.md): ordered implementation plan for producing
  `ariadne-cps/foundation`.

## Historical baseline

The following files describe the baseline captured on 2026-09-29 at
`949d7e044ae65837fc02e6387701b10e7c15ddc6`. They are intentionally retained
as historical evidence and are **not** the current dependency graph:

- [Coupling analysis](coupling-analysis.md)
- [Dependency graph](dependency-graph.svg)
- [Dependency summary](dependency-summary.csv)
- [Include evidence](include-evidence.csv)

At that baseline all ten source directories formed one strongly connected
component with 55 direct dependency relations.

## Current direction

The decoupling work has moved to `main`; the working branch for the next step
is `decouple-foundation`.

Interval has already been removed from `source/geometry` and is consumed from
the standalone `ariadne-cps/interval` repository. The next package boundary is
Foundation: it will contain the current `algebra` and `function` modules and
will depend on `configuration` and `interval`.

The immediate task is therefore not repository movement. It is to remove every
remaining semantic dependency from algebra/function to Ariadne modules that
will stay outside Foundation. In particular, logging is not an intended
Foundation dependency; the one direct include currently found in
`function/calculus_base.hpp` is unused and should be treated as cleanup.
