# Current decoupling state

Status date: 2026-10-03  
Working branch: `decouple-foundation`

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

## Algebra

Algebra has a real dependency on Interval. Current examples include
`Differential<Float*UpperInterval>`, interval-valued expansions, matrices and
sweepers.

The unwanted edge is `algebra -> function`:

- `algebra/graded.hpp` includes `function/procedure.hpp` for
  `compute_procedure`;
- `algebra/algebra_operations.tpl.hpp` includes
  `function/taylor_series.hpp` for TaylorSeries/AnalyticFunction composition;
- `algebra/dense_differential.cpp` includes `function/functional.hpp`, though
  it is not currently part of the algebra object target.

These integrations must move upward so Algebra can be extracted as a repository
built above Interval.

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
