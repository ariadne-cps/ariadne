add_library(ariadne-solvers OBJECT
    solver.cpp
    integrator.cpp
    bounder.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    linear_programming.cpp
    nonlinear_programming.cpp
    runge_kutta_integrator.cpp
)

cotire(ariadne-solvers CONFIGURATIONS None)

if(COVERAGE)
    target_compile_options(ariadne-solvers PUBLIC ${COVERAGE_COMPILER_FLAGS})
endif()