# Migration of the former dynamics directory

The obsolete `source/dynamics` directory has been removed. All internal consumers
now include the owning system, evolution, or reasoning headers directly.

| Ownership | Files |
|---|---|
| `system/` | vector fields, iterated maps, differential inclusions, hybrid automata, mathematical first-order PDE definition and differential-inclusion gridding |
| `evolution/` | enclosures and flow tubes, orbits, configuration interface for evolvers, scaling policy, continuous evolver and simulator implementations, first-order PDE evolution |
| `reasoning/` | reachability analyser, reachability interfaces and template implementation, storage and enclosure-to-storage adapters |
| `dynamics/` | removed; former forwarding headers no longer supported |

The build uses sources from their owning directories. Includes of the historical
`dynamics/` paths are no longer supported; use the owning module's headers.
This is an intentional source-compatibility break for external C++ consumers.

The former `tests/dynamics/` test sources have been reorganized by responsibility:
propagation, enclosure and simulator tests live under `tests/evolution/`, while
finite/infinite-time reachability and safety verification tests live under
`tests/reasoning/`. All twelve original test sources, executable names and CTest
names are preserved; the Windows timeout for the infinite-time reachability test
is also preserved. The PDE regression lives under `tests/evolution/`.

Only `tests/evolution/` and `tests/reasoning/` remain in this repository.
The old model-only `test_differential_inclusion`, `test_iterated_map`, and
`test_vector_field` have been removed here: the `system` repository is the
sole owner of those tests. Their exact historical cases have **not** been
ported to that repository by this change, so some scenario-specific coverage
is intentionally lost in Ariadne; the `system` submodule is unchanged.

**The boundary is about responsibility, not the input type**: numerical
evolution of a PDE remains evolution regardless of whether the input model
contains matrices, symbolic functions or PDE boundary coefficients. A
reachability analyser may depend on an evolver; system definitions must never
include an analyser or an evolver.

See [PDE architecture](system-evolution-pde.md) for the detailed mathematical
model, solver interface, and the deliberate removal of obsolete 1D/2D wave prototypes.
