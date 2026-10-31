# System definitions and PDE evolution

## Dependency boundary

A dynamical **system** specifies a state space and mathematical law. It does
not choose time steps, grid resolutions, integrators, result enclosures or
reachability algorithms. **Evolution** computes solutions or approximations;
**analysis** consumes evolvers and results to answer reachability and safety
questions.

```text
analysis -> evolution -> system -> hybrid -> labelled
```

The `ariadne-system` target must not depend on `evolution`, `analysis` or
the compatibility layer `dynamics/`.

## First-order PDE

[`source/system/first_order_pde.hpp`](../../source/system/first_order_pde.hpp)
defines the mathematical model `FirstOrderPDE`:

```text
A du/dt + sum_j B_j du/dx_j = f(u)
```

It owns `A` and `Bs` and borrows `Ds`, `Ts` and `f` by reference.
Those referenced objects must outlive the model and any solver invocation.
An initial condition, numerical grid, time step and computed solution belong
to evolution, not to the model.

[`source/evolution/first_order_pde.hpp`](../../source/evolution/first_order_pde.hpp)
provides `FirstOrderPDESolution<PR>` (spatial step `h`, temporal step
`tau`, sampled trajectory `uts`, error estimate) and the solver
`first_order_pde`. The solver is implemented in
`source/evolution/first_order_pde.cpp` and consumes the model, initial
condition `phi0` and numeric precision `PR`.

New code should include the layers it actually needs:

```cpp
#include "system/first_order_pde.hpp"      // PDE specification
#include "evolution/first_order_pde.hpp"   // PDE evolution and solution
```

The retained regression is
[`tests/evolution/test_first_order_pde.cpp`](../../tests/evolution/test_first_order_pde.cpp)
and its CTest name is `test_first_order_pde`.

## Removed historical 1D/2D wave prototypes

The historical `1D_pde` and `2D_pde` prototypes have been intentionally
removed, including their model/parameter classes, finite-difference solvers,
associated tests and `dynamics/` compatibility forwarding headers. They
were special-purpose examples, not validated general PDE interfaces; their
discretization assumptions, hardcoded settings, and numerical tests did not
justify keeping them in the core library.

In particular, `WaveEquation1D`, `Parameter1D`, `Parameter2D`,
`pde_1d_solver`, `pde_2d_solver`, `linspace1d`, `linspace2d`,
`set_ic`, `set_ic_2d` and `gaussian_function` from the prototypes are
no longer provided. Compatibility with these APIs was deliberately broken;
the original implementations remain recoverable from Git history.

Retaining `FirstOrderPDE` does not imply a complete generic PDE framework
or new stability and accuracy guarantees. These require independent
mathematical verification and focused numerical tests.
