message(${PROJECT_SOURCE_DIR})

add_library(ariadne-solvers OBJECT
    solver.cc
    integrator.cc
    constraint_solver.cc
    simplex_algorithm.cc
    linear_programming.cc
    nonlinear_programming.cc
    runge_kutta_integrator.cc
)

cotire(ariadne-solvers)
