# Ariadne Foundation decoupling work plan

Status: active  
Branch: `decouple-foundation`  
Historical baseline: `949d7e044ae65837fc02e6387701b10e7c15ddc6` (2026-09-29)  
Current reference: `main` at `7be5bad2d8131ce5e7fec422fe86435a64c69084`  
Last updated: 2026-10-03

See [current-state.md](current-state.md) for the verified current structure and
the evidence behind this plan.

## Goal

Create `ariadne-cps/foundation` from the current algebra and function modules.

The intended external Foundation dependencies are:

1. `configuration`;
2. `interval`.

Foundation must not require Ariadne's remaining geometry, symbolic, io, solving,
dynamics, hybrid, logging or threading modules.

Repository extraction is the final step. The boundary must first be enforced
inside Ariadne.

## Current facts

- Interval is already a standalone repository and has been removed from
  `source/geometry`.
- Current algebra interval usage goes directly through `interval/...`; the old
  baseline `algebra -> geometry` interval edge is obsolete.
- Algebra's remaining upward coupling is concentrated in function-specific
  integrations.
- Function still has real geometry and symbolic dependencies.
- No semantic algebra/function dependency on logging has been established.
  `function/calculus_base.hpp` contains one unused logging include.
- The current top-level Ariadne build still aggregates logging/threading objects;
  aggregate linkage is not evidence that Foundation itself needs those
  dependencies.

## Milestones

### F0 - Make the documentation describe the current graph

- [x] Record the standalone Interval boundary.
- [x] Separate the 2026-09-29 baseline from the current repository state.
- [x] Record Foundation as `algebra + function` with external dependencies
      `configuration + interval`.
- [x] Record the currently observed residual algebra/function edges.
- [ ] Add a repeatable dependency scanner so the current graph is generated
      rather than maintained manually.

Exit criterion: contributors can distinguish historical coupling evidence from
the active Foundation plan.

### F1 - Make algebra independent of function

- [ ] Move `compute_procedure` and its `Procedure` dependency out of
      `algebra/graded.hpp`.
- [ ] Move TaylorSeries/AnalyticFunction composition code out of
      `algebra/algebra_operations.tpl.hpp`.
- [ ] Decide whether to delete, revive, or relocate
      `algebra/dense_differential.cpp`; it is not in the current algebra
      object target but still carries a function include.
- [ ] Compile every public algebra header without adding function include paths.
- [ ] Give algebra an explicit target-level dependency on its actual lower
      dependencies instead of inheriting Ariadne-wide linkage.

Exit criterion: no algebra source/header includes `function/` and algebra
builds/tests without the function target.

### F2 - Remove accidental infrastructure dependencies from function

- [ ] Remove the unused `logging/logging.hpp` include from
      `function/calculus_base.hpp`.
- [ ] Verify no logging symbol is required by any function implementation.
- [ ] Verify function and algebra do not require threading directly.
- [ ] Stop relying on parent-directory `link_libraries(threading interval)` for
      the future Foundation targets; declare only real target dependencies.

Exit criterion: logging and threading can be absent from a Foundation-only
configure/build.

### F3 - Establish a primitive domain boundary

- [ ] Design a low-level Box/domain value type that depends only on the lower
      Foundation/Interval stack.
- [ ] Split that primitive from the current `geometry/box.hpp`, whose public
      surface currently also depends on Point and SetInterface.
- [ ] Move function domain declarations to the primitive Box contract.
- [ ] Leave Box set algorithms, function integrations and drawing adapters above
      the primitive boundary.
- [ ] Re-run public-header compilation for function after removing
      `geometry/box*.hpp` dependencies.

Exit criterion: the core Function API can represent scalar/vector domains
without including Ariadne geometry.

### F4 - Move geometry integrations above function core

- [ ] Move measurable-function integrations that depend on
      `geometry/set.hpp`, `measurable_set.hpp` and `set_wrapper.hpp` out of
      the Foundation function core.
- [ ] Move multifunction/Taylor-multifunction integrations that depend on
      FunctionSet/SetWrapper out of the core, retaining only genuinely generic
      function/multifunction abstractions.
- [ ] Decide whether these adapters belong to geometry or a separate integration
      layer.

Exit criterion: no Foundation function source/header includes `geometry/`.

### F5 - Split low-level symbolic machinery from symbolic integration

- [ ] Move concrete Expression/Function conversion code out of
      `function/function.cpp` to the symbolic side.
- [ ] Classify `symbolic/templates.hpp`: it currently depends only on numeric
      and is used by Formula/Procedure, so evaluate moving it to Foundation under
      a non-symbolic-specific name/location.
- [ ] Classify `symbolic/constant.hpp` and `symbolic/identifier.hpp` as either
      low-level Foundation value types or symbolic-owned types; avoid keeping a
      path-level dependency merely for historical naming.
- [ ] Ensure Foundation has no dependency on the high-level Expression, Space or
      Variable implementation.

Exit criterion: no Foundation target includes high-level `symbolic/` headers;
any retained low-level machinery is owned by Foundation itself.

### F6 - Standalone Foundation target

- [ ] Create component-level CMake targets for algebra and function with explicit
      public/private dependencies.
- [ ] Add public-header compilation checks.
- [ ] Make Foundation tests link only to Foundation targets.
- [ ] Add standalone Unix/Windows/Coverage CI.
- [ ] Add an installed/external-consumer test.
- [ ] Confirm Foundation configures without Ariadne source directories.

Exit criterion: Foundation is independently buildable and testable with only
its declared external dependencies.

### F7 - Repository extraction

- [ ] Create/populate `ariadne-cps/foundation`.
- [ ] Depend on the required `configuration` and `interval` revisions.
- [ ] Move algebra/function sources, tests and package configuration.
- [ ] Replace Ariadne's in-tree algebra/function copies with the Foundation
      dependency.
- [ ] Keep Ariadne integration CI spanning Foundation plus the remaining
      high-level modules.

Exit criterion: Ariadne consumes Foundation as a versioned external repository
and no duplicate algebra/function implementation remains in Ariadne.

## Immediate implementation slice

The smallest useful code slice after this documentation update is:

1. remove the unused logging include;
2. move the two function-specific public-template integrations out of algebra;
3. add an algebra public-header compile check that runs without function;
4. then tackle Box as the first structural function/geometry boundary.

This order deliberately proves the easy boundary first. Starting with Box would
mix a local cleanup problem with the larger geometry split and make failures
harder to attribute.

## Decisions

| ID | Date | Decision | Rationale |
|---|---|---|---|
| D-001 | 2026-09-29 | Preserve the original coupling analysis as a historical baseline. | It remains useful evidence but no longer describes the repository after extraction work. |
| D-002 | 2026-10-03 | Build the next Foundation package from algebra and function, depending externally on configuration and interval. | Interval is already standalone; algebra/function are the next low-level reusable layer. |
| D-003 | 2026-10-03 | Do not treat logging as a Foundation dependency. | The only direct function include found is unused; aggregate Ariadne linkage must not define the package boundary. |
| D-004 | 2026-10-03 | Treat Box as a primitive-boundary refactor, not a wholesale geometry move. | Current Box headers/implementation mix primitive domain representation with Point/SetInterface/function/io responsibilities. |
| D-005 | 2026-10-03 | Separate low-level symbolic templates from high-level Expression/Function bridges. | The former may belong below symbolic; the latter are integration code and should stay above Foundation. |

## Validation discipline

For each boundary change:

1. compile the changed public headers in isolation;
2. build the affected component target without relying on transitive sibling
   include directories;
3. run its focused tests;
4. run the Ariadne integration build/tests;
5. record the removed/added dependency edge in this document or regenerated
   dependency data.

Do not interpret successful compilation through an aggregate Ariadne target as
proof that a standalone package boundary is correct.
