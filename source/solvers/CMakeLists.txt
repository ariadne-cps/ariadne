add_library(ariadne-solvers OBJECT
    solver.cpp
    integrator.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    linear_programming.cpp
    nonlinear_programming.cpp
    runge_kutta_integrator.cpp
)

cotire(ariadne-solvers)
