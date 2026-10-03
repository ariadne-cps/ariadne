# Current decoupling state

Status date: 2026-10-03  
Ariadne reference: `main` at `7be5bad2d8131ce5e7fec422fe86435a64c69084`  
Working branch: `decouple-foundation`

This document describes the current repository/package structure. The dependency
CSV, SVG and include evidence in this directory remain the historical baseline
captured on 2026-09-29 and must not be read as the current graph.

## Current repository structure

Ariadne no longer owns the interval implementation under `source/geometry`.
It consumes `ariadne-cps/interval` as a pinned submodule.

The relevant standalone dependency chain is currently:

```text
configuration

interval
  -> numeric
       -> paradigm
            -> utility
       -> configuration
```

`paradigm` also carries its own configuration/utility support and Python
support as repository dependencies. Ariadne separately carries `threading`,
but that is an Ariadne-level integration dependency, not part of the intended
Foundation contract.

The source directories still owned by Ariadne are:

```text
algebra
function
solving
geometry
dynamics
symbolic
io
hybrid
```

## Target Foundation boundary

The next extraction target is `ariadne-cps/foundation`.

Foundation is intended to contain the current `algebra` and `function`
modules and to depend externally on:

1. `configuration`;
2. `interval`.

The lower numeric/paradigm/utility chain is reached through `interval`; it is
not a reason to keep copies of those modules in Ariadne.

Before extraction, `algebra` and `function` must have no semantic dependency
on modules that remain in Ariadne.

## Current residual dependencies

### Algebra

The old baseline reported an `algebra -> geometry` dependency because Interval
was still stored under geometry. That dependency is gone in the current source:
the relevant algebra files now include `interval/interval.hpp` or
`interval/interval.decl.hpp` directly.

No direct algebra dependency on logging, threading, io, solving, dynamics,
hybrid or symbolic was found in the current source inspected for this step.

The remaining upward dependency is `algebra -> function`:

- `algebra/graded.hpp` includes `function/procedure.hpp` for
  `compute_procedure`;
- `algebra/algebra_operations.tpl.hpp` includes
  `function/taylor_series.hpp` and implements composition involving
  `TaylorSeries` / `AnalyticFunction`;
- `algebra/dense_differential.cpp` includes `function/functional.hpp`, but
  this file is not part of the current algebra CMake object target.

These integrations should move to the function side (or to a Foundation
integration header owned above the algebra core). The algebra core should not
know function types.

### Function -> logging

There is no demonstrated semantic logging dependency.

The only current direct logging include found in algebra/function is:

```text
source/function/calculus_base.hpp -> logging/logging.hpp
```

No logging symbol is used in that header. Treat this as an obsolete include to
remove, not as a Foundation dependency. The aggregate Ariadne libraries still
contain logging objects through other top-level dependencies; that must not be
mistaken for a module-level dependency of Foundation.

### Function -> geometry

This is the largest remaining boundary problem.

Current direct examples are:

- `function/domain.hpp -> geometry/box.hpp`;
- `function/function.decl.hpp -> geometry/box.decl.hpp`;
- `function/measurable_function.hpp -> geometry/set.hpp`,
  `geometry/measurable_set.hpp`, `geometry/set_wrapper.hpp`;
- `function/multifunction.hpp -> geometry/set.hpp`, `geometry/box.hpp`,
  `geometry/function_set.hpp`, `geometry/set_wrapper.hpp`;
- `function/taylor_multifunction.* -> geometry/box.hpp` and
  `geometry/set_wrapper.hpp`.

The `Box` case is architectural, not just an include cleanup. The current
`geometry/box.hpp` still depends on `geometry/point.hpp` and
`geometry/set_interface.hpp`, while `geometry/box.cpp` mixes function,
algebra and io integrations. Therefore moving `geometry/box.*` wholesale into
Foundation would import unwanted geometry responsibilities.

The likely direction is to split a low-level Box/domain value type from
higher-level set, algorithm and drawing integrations. That primitive can then
sit with Interval/Foundation, while geometry retains the higher-level behaviour.

The measurable-function and multifunction set integrations should instead move
upward out of the Foundation function core unless a smaller generic contract is
identified.

### Function -> symbolic

This edge also contains two different concerns and should not be handled as one
bulk move.

Low-level candidates:

- `symbolic/templates.hpp` depends only on numeric headers in the current
  source and is used by `function/formula.hpp` and
  `function/procedure.hpp`. It is a candidate for relocation to a lower,
  non-symbolic-specific support layer inside Foundation.
- `symbolic/constant.hpp` additionally depends on
  `symbolic/identifier.hpp`; both are small, but their ownership should be
  decided from semantics rather than path names.

High-level bridges:

- `function/function.cpp` includes `symbolic/expression.hpp` specifically
  for conversions between Expression and Function.
- Expression/space/variable conversion code should remain above Foundation,
  most naturally on the symbolic side, rather than forcing Foundation to depend
  on the full symbolic module.

## Proposed decoupling order

1. Remove the unused logging include from `function/calculus_base.hpp` and
   ensure Foundation does not link logging merely because Ariadne does.
2. Remove `algebra -> function` by moving Procedure/TaylorSeries-specific
   integrations out of the algebra core.
3. Split the primitive Box/domain representation from the current geometry
   implementation and make function depend only on that primitive.
4. Move measurable/multifunction set integrations out of the function core.
5. Move concrete Expression/Function conversion code to the symbolic side.
6. Relocate or rename genuinely low-level symbolic templates/constants that are
   needed by Formula/Procedure.
7. Give algebra and function explicit standalone CMake targets whose only
   external package dependencies are the intended Foundation dependencies.
8. Add public-header compilation tests and a minimal external consumer for the
   future Foundation target.
9. Only then move algebra/function into `ariadne-cps/foundation` and replace
   their Ariadne copies with the Foundation dependency.

## Validation criterion

The boundary is ready when a standalone Foundation checkout can configure,
compile public headers, build its implementation, run its tests, and satisfy a
small external consumer with only its declared dependencies, without any
include or link path into Ariadne's `geometry`, `symbolic`, `io`,
`solving`, `dynamics`, `hybrid`, `logging` or `threading` modules.
