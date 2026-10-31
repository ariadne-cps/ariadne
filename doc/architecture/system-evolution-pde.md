# System definitions and PDE evolution

## Dependency boundary

A dynamical **system** specifies a state space and a mathematical law. It must
not choose a time step, a grid resolution, an integrator, an enclosure, or a
reachability algorithm. **Evolution** consumes that law and initial data to
compute trajectories or approximations; **analysis** consumes evolvers to
answer reachability and safety questions.

```text
analysis -> evolution -> system -> hybrid -> labelled
```

The CMake sources of `ariadne-system` must not include `evolution/`,
`analysis/`, or `dynamics/`. Compatibility headers in `dynamics/` are for
old clients only and must not be used by new internal code.

## First-order PDE

**Model**: `source/system/first_order_pde.hpp` declares `FirstOrderPDE`,
the system

```text
A du/dt + sum_j B_j du/dx_j = f(u)
```

with coefficient matrices `A`, `Bs`, diagonalizations `Ds`, `Ts`, and
reaction/source function `f`. The `Ds`, `Ts`, and `f` members are borrowed
references: all three referenced values must outlive the model and any solver
call that uses it. This is the **existing ownership contract**; no
discretization or result data is stored in the model.

**Evolution**: `source/evolution/first_order_pde.hpp` declares
`FirstOrderPDESolution<PR>` (spatial step `h`, temporal step `tau`,
sampled trajectory `uts`, and error bound) and the `first_order_pde`
algorithm. Its implementation is in `source/evolution/first_order_pde.cpp`.
The initial condition `phi0` and precision `PR` are solver inputs, not
properties of the PDE. Numerical interpolation, step-count checks, and error
estimation belong entirely to evolution.

In particular, `first_order_pde` **computes an approximation/solution**;
it does not perform reachability analysis. Keeping it in `analysis` would
reverse the dependency boundary.

## Legacy one- and two-dimensional wave prototypes

The original `1D_pde.hpp` and `2D_pde.hpp` combined physical data and a
finite-difference algorithm in one header. They are now split as follows:

| Definition (system) | Sampling and propagation (evolution) |
| --- | --- |
| `system/one_dimensional_wave_equation.hpp`: `WaveEquation1D<PR>` | `evolution/one_dimensional_wave_solver.hpp`: `Parameter1D<PR>`, `set_ic`, `linspace1d`, `pde_1d_solver` |
| `system/two_dimensional_wave_equation.hpp`: `Parameter2D<PR>` | `evolution/two_dimensional_wave_solver.hpp`: `set_ic_2d`, `linspace2d`, `pde_2d_solver`, `gaussian_function` |

`Parameter1D` remains source compatible with old callers and extends the
physical `WaveEquation1D` model with `CourantNumber`. This latter quantity
is **a finite-difference tuning parameter**, not part of the physical
definition. Keeping it only in the evolution-side compatibility configuration
prevents it from contaminating the new system definition.

The 2D `Parameter2D` represents physical information for one spatial axis.
The existing 2D routine still chooses its propagation speed, final time and
time step internally. Both historical solvers also choose particular initial
profiles and boundary conditions. **They are examples of finite-difference
schemes, not a general or rigorously validated PDE solver interface.**
No claim of a general PDE framework or reliable stability/error guarantees
is implied by moving them. Refactoring their numerical algorithms and
validation must be separate, behavior-tested work.

The unused `1D_pde.cpp` and `2D_pde.cpp` copies were removed; they had
different signatures and were not part of the CMake build. The functions
actually tested were the templates in the headers.

## Compatibility and usage

Old source paths continue to work via forwarding includes:

```cpp
#include "dynamics/first_order_pde.hpp"
#include "dynamics/1D_pde.hpp"
#include "dynamics/2D_pde.hpp"
```

New code should express intent explicitly:

```cpp
#include "system/first_order_pde.hpp"            // mathematical PDE
#include "evolution/first_order_pde.hpp"         // solution type + algorithm
#include "system/one_dimensional_wave_equation.hpp"
#include "evolution/one_dimensional_wave_solver.hpp"
```

The PDE regression is `tests/evolution/test_pdes.cpp` and runs as
`ctest -R '^test_pdes$'`. It tests existing PDE behavior; relocating the
interfaces does **not** constitute new solver accuracy or stability testing.
