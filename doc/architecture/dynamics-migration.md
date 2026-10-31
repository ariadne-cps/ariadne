# Migration of the former dynamics directory

The obsolete `source/dynamics` directory is no longer a build unit. Its
remaining headers are forwarding adapters for source compatibility.

| Ownership | Files |
|---|---|
| `system/` | vector fields, iterated maps, differential inclusions, hybrid automata, mathematical first-order PDE definition |
| `evolution/` | enclosures and flow tubes, orbits, configuration interface for evolvers, scaling policy, continuous evolver and simulator implementations, differential inclusion gridding, first-order PDE evolution |
| `analysis/` | reachability analyser, reachability interfaces and template implementation, storage used in reachability analysis |
| `dynamics/` | forwarding headers only, to be deprecated in a future API-breaking release |

The build uses sources from their owning directories. A compatibility
`#include "dynamics/name.hpp"` is permitted for existing consumers but no
new lower-level implementation should depend on these forwarding headers.

The legacy `tests/dynamics/` executable organization is retained for now to
preserve the original behavior-focused regression and its CTest names; the
PDE regression itself now lives under `tests/evolution/test_first_order_pde.cpp`.
A subsequent change can relocate other tests without conflating test
ownership with production dependency direction.

**The boundary is about responsibility, not the input type**: numerical
evolution of a PDE remains evolution regardless of whether the input model
contains matrices, symbolic functions or PDE boundary coefficients. A
reachability analyser may depend on an evolver; system definitions must never
include an analyser or an evolver.

See [PDE architecture](system-evolution-pde.md) for the detailed mathematical
model, solver interface, and the deliberate removal of obsolete 1D/2D wave prototypes.
