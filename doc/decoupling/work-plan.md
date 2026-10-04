# Ariadne decoupling work plan

Status: active  
Branch: `decouple-algebra`  
Historical baseline: `949d7e044ae65837fc02e6387701b10e7c15ddc6`  
Last updated: 2026-10-04

## Target chain

```text
utility -> foundation -> numeric -> interval -> algebra -> function
```

Foundation is the renamed former Paradigm package, not an Algebra/Function
aggregator.

## F0 - Rename Paradigm package to Foundation

- [x] Rename GitHub repository `ariadne-cps/paradigm` to
      `ariadne-cps/foundation`.
- [x] Rename package/build identifiers: project, target, installed library,
      `PARADIGM_SRC -> FOUNDATION_SRC`.
- [x] Move public headers from `include/paradigm/` to
      `include/foundation/`.
- [x] Preserve `paradigm.hpp` and C++ `Paradigm*` concepts.
- [x] Update Numeric to consume `submodules/foundation`, target
      `foundation`, `FOUNDATION_SRC`, Python Foundation targets, and
      `foundation/...` includes.
- [x] Layer the Python binding interfaces as
      `pyariadne-foundation -> pyariadne-numeric -> pyariadne-interval -> pyariadne-algebra -> pyariadne`.
- [x] Make each Python layer consume only the aggregate binding
      interface/header of its direct lower-level repository, with no manual
      references to nested submodule binding paths.
- [x] Start reducing Ariadne `python/bindings/utilities.hpp` to
      Ariadne-specific helpers by removing copies already owned by
      python-common/Foundation/Numeric/Interval.
- [ ] Promote the coordinated Foundation/Numeric revisions and cascade updated
      pins through Interval and Ariadne.

## F1 - Make Algebra independent of Function

- [x] Move `compute_procedure` out of `algebra/graded.hpp`.
- [x] Remove TaylorSeries-specific composition from
      `algebra/algebra_operations.tpl.hpp`; retain `AnalyticFunction`, which is
      defined by Algebra itself.
- [x] Remove the unnecessary Function include from
      `algebra/dense_differential.cpp`.
- [x] Compile Algebra public headers without Function as part of a completely
      green standalone validation.
- [x] Give Algebra an explicit target-level dependency on Interval.
- [x] Create the standalone `ariadne-cps/algebra` repository above Interval.
- [x] Move `Polynomial`, `UnivariateChebyshevPolynomial` and
      `MultivariateChebyshevPolynomial` from Function into Algebra, including
      focused tests and Python bindings.
- [x] Complete the Chebyshev algebra-operation contract with unary `Pos`.
- [x] Move `SweeperBase` and `RelativeSweeperBase` implementation code out
      of `function/taylor_model.tpl.hpp` and into the Algebra repository.
- [x] Get all standalone Algebra CI jobs green: C++, Python, installation and
      external-consumer/tutorial checks.
- [x] Update Ariadne to consume the validated standalone Algebra revision and
      remove the local Algebra implementation, duplicate bindings, and Algebra-owned tests/demonstrations.

Exit criterion: Algebra depends only on its intended lower stack, with Interval
as its immediate repository dependency, and the standalone Algebra CI is green.

## F2 - Remove accidental Function infrastructure dependencies

- [ ] Remove unused `logging/logging.hpp` from
      `function/calculus_base.hpp`.
- [ ] Verify Function has no direct Logging or Threading requirement.
- [ ] Use explicit target dependencies rather than Ariadne-wide linkage.

## F3 - Remove Function -> Geometry

- [ ] Split a low-level Box/domain value type from current geometry Box.
- [ ] Move Function domain declarations to that primitive.
- [ ] Move measurable-function and multifunction set integrations upward.
- [ ] Ensure Function public headers no longer include `geometry/`.

## F4 - Remove Function -> high-level Symbolic

- [ ] Move Expression/Function conversion code to the Symbolic side.
- [ ] Classify low-level symbolic templates/constants by semantics and relocate
      only those genuinely required below Symbolic.
- [ ] Ensure Function no longer depends on high-level Expression/Space/Variable
      implementation.

## F5 - Extract Function

- [ ] Add standalone Function CMake target and public-header checks.
- [ ] Add focused tests and external-consumer test.
- [ ] Extract `ariadne-cps/function` above Algebra.
- [ ] Replace Ariadne in-tree copies with pinned repository dependencies.

## Validation discipline

For every boundary change:

1. compile affected public headers in isolation;
2. build the component without sibling include leakage;
3. for Python bindings, consume only the aggregate `pyariadne-<direct-dependency>`
   interface and direct dependency's aggregate Python header; never add include
   paths or source references to nested submodules;
4. keep each layer's Python utility header minimal: move or remove helpers that
   are already owned by a lower layer instead of duplicating them;
5. run focused tests;
6. run Ariadne integration tests;
7. record the removed dependency edge here.

Do not treat successful aggregate Ariadne linkage as proof of a valid standalone
boundary. In particular, a Python build is not considered decoupled if it works
only because transitive repositories have been exposed manually through nested
submodule paths.
