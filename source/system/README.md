# System and evolution layering

System code contains definitions only: maps, vector fields, differential
inclusions, hybrid automata, and PDE problem specifications. No propagation,
numerical discretization, trajectory computation, or reachability belongs
here. All system implementations must remain independent of the higher
`evolution` and `analysis` layers.

See [System definitions and PDE evolution](../../doc/architecture/system-evolution-pde.md)
for the boundary, PDE model ownership, compatibility paths, and caveats.
