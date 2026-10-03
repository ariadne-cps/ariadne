# Decoupling Ariadne components

This directory records the architectural analysis and current work for moving
low-level Ariadne components into independently buildable repositories.

Dependency arrows in the historical analysis use **provider -> dependent**.

## Read this first

- [Current state](current-state.md): verified repository structure as of 2026-10-03
  and the active extraction chain.
- [Working plan](work-plan.md): ordered work for Algebra and Function decoupling.

## Historical baseline

The following files describe the baseline captured on 2026-09-29 at
`949d7e044ae65837fc02e6387701b10e7c15ddc6`. They are retained as historical
evidence and are not the current graph:

- [Coupling analysis](coupling-analysis.md)
- [Dependency graph](dependency-graph.svg)
- [Dependency summary](dependency-summary.csv)
- [Include evidence](include-evidence.csv)

## Current direction

The old `ariadne-cps/paradigm` repository has been renamed to
`ariadne-cps/foundation`. Foundation remains the low-level logical/computational
paradigm package; it is not an aggregator for Algebra and Function.

The intended repository chain is now:

```text
utility
  -> foundation
       -> numeric
            -> interval
                 -> algebra
                      -> function
```

Configuration remains shared build infrastructure where required.

The same one-level dependency rule applies to the Python bindings. Each repository
exports an aggregate `pyariadne-<component>` interface and its component-level
Python utility/header surface. A higher layer consumes only the Python interface
of its **direct** lower-level repository; it must not add include paths or source
references into nested submodules. Lower-level Python binding requirements are
propagated transitively by the direct dependency.

As part of this work, Ariadne's `python/bindings/utilities.hpp` is being reduced
to Ariadne-specific helpers only. Helpers already owned by python-common,
Foundation, Numeric or Interval are included through the direct Interval Python
surface rather than copied or redefined in Ariadne.

Interval is already standalone. The next architectural task is to make Algebra
independent of Function, extract Algebra above Interval, then make Function
independent of the Ariadne modules that must remain above it before extracting
Function above Algebra.

The direct `logging/logging.hpp` include in
`function/calculus_base.hpp` is currently unused and is cleanup, not a package
dependency.
