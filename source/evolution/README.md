# Evolution layer

Evolution contains simulators, evolvers, enclosures/orbits, and numerical
propagation schemes, including first-order PDE evolution. It depends on
`system` but not on `analysis`. See
[System definitions and PDE evolution](../../doc/architecture/system-evolution-pde.md)
for the model/solver separation and the deliberate removal of obsolete wave prototypes.

The grid-based `Storage` and `Orbit<Storage>` implementations reside in
`analysis/`; the evolution object library does not include or link to analysis.
Two legacy `Enclosure` methods still mention `Storage` in their signatures
for compatibility, but are defined in `analysis/enclosure_storage.cpp`.
New evolution APIs must not expose analysis storage types.
The evolution CTest dependency guard checks direct include regressions.
The obsolete `source/dynamics/` forwarding headers have been removed.
