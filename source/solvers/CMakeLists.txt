set(LIBRARY_NAME ariadne-solvers)

add_library(${LIBRARY_NAME} OBJECT
    solver.cpp
    smt_solver.cpp
    integrator.cpp
    inclusion_integrator.cpp
    bounder.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    linear_programming.cpp
    nonlinear_programming.cpp
    runge_kutta_integrator.cpp
)


