# Current decoupling state

Status date: 2026-10-04  
Working branch: `decouple-algebra`

The CSV, SVG and include evidence in this directory remain the historical
2026-09-29 baseline. This file describes the active architecture.

## Repository chain

The low-level repositories are being organised as:

```text
configuration   (shared build support)

utility
  -> foundation
       -> numeric
            -> interval
                 -> algebra
                      -> function
```

`foundation` is the renamed former `paradigm` repository. Its public header
directory becomes `include/foundation`, while `paradigm.hpp` and the C++
concepts `Paradigm`, `ParadigmCode`, `ParadigmTraits`, etc. retain their
semantic names.

Numeric consumes Foundation as `submodules/foundation`, links the
`foundation` target, and aggregates `FOUNDATION_SRC`.

## Python binding layering

Python bindings follow the same repository layering as the C++ libraries:

```text
pyariadne-foundation
  -> pyariadne-numeric
       -> pyariadne-interval
            -> pyariadne-algebra
                 -> pyariadne
```

Each layer has two responsibilities:

1. own and compile the Python bindings for that component;
2. expose an aggregate Python interface for the next layer, including the
   binding include requirements inherited from its direct lower dependency.

A consumer must know only its **direct** lower-level Python dependency. It must
not reach through that dependency to nested repositories with paths such as
`submodules/.../submodules/...`, nor manually add binding include directories
for transitive components.

Concretely, Ariadne consumes the Interval Python binding surface only. Its local
`python/bindings/utilities.hpp` includes the aggregate Interval utility header,
which in turn obtains Numeric and Foundation Python support through the lower
layers. Ariadne must therefore contain no direct Python include-path knowledge of
Numeric or Foundation.

The same rule applies recursively: Interval consumes Numeric's Python surface;
Numeric consumes Foundation's Python surface. The public `pyariadne-<component>`
INTERFACE target is responsible for propagating what the next layer needs.

This layering is also being used to shrink
`ariadne/python/bindings/utilities.hpp`. It must contain only helpers belonging
to the Ariadne layer. Generic Python machinery belongs in python-common;
representation helpers belong in Foundation; numeric operators and arithmetic
binding helpers belong in Numeric; interval-specific representations belong in
Interval. Duplicating these definitions in Ariadne is a boundary violation and
can also produce C++ redefinition errors once the lower aggregate headers are
correctly visible.

## Ariadne aggregate targets

The former `ariadne-core` and `ariadne-kernel` shared-library aggregates have
been removed on `decouple-algebra`. In-tree tests and benchmarks now link the
single `ariadne` library directly. The corresponding Python binding object
split has also been removed: Ariadne-owned bindings are compiled through one
`pyariadne-bindings-obj` target, while `pyariadne-module-obj` remains separate
only for the Python module entry point.

This deliberately favours a single aggregate at the Ariadne level while the
repository stack is being split into independently packaged lower layers.

## Algebra

Algebra has a real dependency on Interval. Current examples include
`Differential<Float*UpperInterval>`, interval-valued expansions, matrices and
sweepers.

The former unwanted edge `algebra -> function` has been removed on
`decouple-algebra`:

- `compute_procedure` was removed from `algebra/graded.hpp` and placed next
  to its only consumer in the solving integrator;
- TaylorSeries-specific composition and the
  `function/taylor_series.hpp` include were removed from
  `algebra/algebra_operations.tpl.hpp`; the remaining
  `AnalyticFunction` composition is algebraic because `AnalyticFunction` is
  defined in `algebra/series.hpp`;
- the unnecessary `function/functional.hpp` include was removed from
  `algebra/dense_differential.cpp`;
- the `ariadne-algebra` target now declares Interval explicitly.

Algebra source files therefore have no direct Function include.

The standalone repository `ariadne-cps/algebra` now exists above Interval on
the coordinated `decouple-algebra` branch. The extraction also moved the
polynomial representations that semantically belong to Algebra:

- `Polynomial` and its implementation/templates were moved from
  `source/function` into Algebra;
- `UnivariateChebyshevPolynomial` and
  `MultivariateChebyshevPolynomial` were moved with their C++ tests and
  Python bindings;
- the Chebyshev operation set was completed with unary `Pos`, matching the
  `DispatchAlgebraOperations` contract already satisfied by `Polynomial`;
- the `SweeperBase` and `RelativeSweeperBase` implementations, which were
  still physically defined in `function/taylor_model.tpl.hpp`, were moved
  into Algebra so that sweeper vtables no longer require Function-owned
  implementation code.

Standalone Algebra CI is green and the extraction has been merged to
`ariadne-cps/algebra:main`. Ariadne on `decouple-algebra` now consumes Algebra
as its direct repository dependency, with Interval supplied through Algebra.
The former local `source/algebra` tree, Function-owned Polynomial and
ChebyshevPolynomial implementations, Sweeper implementations, duplicate Algebra
Python bindings, and Algebra-owned tests/demonstrations have been removed from
Ariadne.

## Function

Function should then be extracted as a repository built above Algebra.

The remaining upward dependencies are mainly:

- Geometry: Box/domain and set/multifunction integrations.
- Symbolic: generic symbolic helpers plus concrete Expression/Function bridges.
- Logging: one unused include in `function/calculus_base.hpp`, not a semantic
  dependency.

The Box case requires a real boundary split because the current geometry Box
mixes primitive domain representation with Point/SetInterface and higher-level
integration code.

## Validation

A repository boundary is ready only when the component configures, builds its
public headers and implementation, runs focused tests, and is consumable
externally using only declared lower-level dependencies.

For Python bindings, validation additionally requires that the component builds
when given only the aggregate Python interface of its direct lower repository.
A successful build that relies on manually exposed nested-submodule include
paths does not validate the boundary.
