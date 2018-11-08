add_library(ariadne-solvers OBJECT
    solver.cpp
    integrator.cpp
    inclusion_integrator.cpp
    bounder.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    linear_programming.cpp
    nonlinear_programming.cpp
    runge_kutta_integrator.cpp
)

cotire(ariadne-solvers CONFIGURATIONS None)

if(COVERAGE)
    target_link_libraries(ariadne-solvers PUBLIC coverage_config)
endif()