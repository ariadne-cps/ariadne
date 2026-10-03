# Ariadne decoupling work plan

Status: active  
Branch: `decouple-foundation`  
Historical baseline: `949d7e044ae65837fc02e6387701b10e7c15ddc6`  
Last updated: 2026-10-03

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
- [ ] Promote the coordinated Foundation/Numeric revisions and cascade updated
      pins through Interval and Ariadne.

## F1 - Make Algebra independent of Function

- [ ] Move `compute_procedure` out of `algebra/graded.hpp`.
- [ ] Move TaylorSeries/AnalyticFunction composition out of
      `algebra/algebra_operations.tpl.hpp`.
- [ ] Resolve `algebra/dense_differential.cpp`.
- [ ] Compile Algebra public headers without Function.
- [ ] Give Algebra explicit target-level dependencies.
- [ ] Extract `ariadne-cps/algebra` above Interval.

Exit criterion: Algebra depends only on its intended lower stack, with Interval
as its immediate repository dependency.

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
3. run focused tests;
4. run Ariadne integration tests;
5. record the removed dependency edge here.

Do not treat successful aggregate Ariadne linkage as proof of a valid standalone
boundary.
