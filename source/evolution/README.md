# Evolution layer

Evolution contains simulators, evolvers, enclosures/orbits, and numerical
propagation schemes, including first-order PDE evolution. It depends on
`system` but not on `reasoning`. See
[System definitions and PDE evolution](../../doc/architecture/system-evolution-pde.md)
for the model/solver separation and the deliberate removal of obsolete wave prototypes.

The grid-based storage and storage orbit implementations reside in `reasoning/`.
The reasoning-only adapter `reasoning/enclosure_storage.hpp` declares the
`outer_approximation(enclosure, grid, fineness)` and
`adjoin_outer_approximation_to(enclosure, storage, fineness)` free functions.
`Enclosure` does not expose reasoning storage types, and the evolution object
library neither includes nor links to `reasoning`.
The evolution CTest dependency guard checks direct include regressions.
The obsolete `source/dynamics/` forwarding headers have been removed.
