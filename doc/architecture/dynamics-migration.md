# Migration of the former dynamics directory

The obsolete `source/dynamics` directory has been removed. All internal consumers
now include the owning system, evolution, or analysis headers directly.

| Ownership | Files |
|---|---|
| `system/` | vector fields, iterated maps, differential inclusions, hybrid automata, mathematical first-order PDE definition and differential-inclusion gridding |
| `evolution/` | enclosures and flow tubes, orbits, configuration interface for evolvers, scaling policy, continuous evolver and simulator implementations, first-order PDE evolution |
| `analysis/` | reachability analyser, reachability interfaces and template implementation, storage used in reachability analysis |
| `dynamics/` | removed; former forwarding headers no longer supported |

The build uses sources from their owning directories. Includes of the historical
`dynamics/` paths are no longer supported; use the owning module's headers.
This is an intentional source-compatibility break for external C++ consumers.

The former `tests/dynamics/` test sources have been reorganized by responsibility:
propagation, enclosure and simulator tests live under `tests/evolution/`, while
finite/infinite-time reachability and safety verification tests live under
`tests/analysis/`. All twelve original test sources, executable names and CTest
names are preserved; the Windows timeout for the infinite-time reachability test
is also preserved. The PDE regression lives under `tests/evolution/`.

Three model-only regression tests (`test_differential_inclusion`,
`test_iterated_map`, `test_vector_field`) temporarily live under
`tests/system_legacy/`. At the currently pinned standalone `system` revision,
equivalent tests for all their construction, missing-dynamics and input
validation cases could not be established. They remain active in Ariadne to
avoid silently dropping coverage, while the external `system` submodule and
its tests remain unchanged.

**The boundary is about responsibility, not the input type**: numerical
evolution of a PDE remains evolution regardless of whether the input model
contains matrices, symbolic functions or PDE boundary coefficients. A
reachability analyser may depend on an evolver; system definitions must never
include an analyser or an evolver.

See [PDE architecture](system-evolution-pde.md) for the detailed mathematical
model, solver interface, and the deliberate removal of obsolete 1D/2D wave prototypes.
