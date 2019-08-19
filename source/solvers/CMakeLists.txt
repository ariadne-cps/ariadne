add_library(ariadne-solvers OBJECT
    solver.cpp
    integrator.cpp
    bounder.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    linear_programming.cpp
    nonlinear_programming.cpp
    runge_kutta_integrator.cpp
    quadratic_programming.cpp
)

cotire(ariadne-solvers CONFIGURATIONS None)

if(COVERAGE)
<<<<<<< HEAD
    target_link_libraries(ariadne-solvers PUBLIC coverage_config)
=======
    target_compile_options(ariadne-solvers PUBLIC ${COVERAGE_COMPILER_FLAGS})
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
endif()
