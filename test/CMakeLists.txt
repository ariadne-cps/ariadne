set(UNIT_TESTS_KERNEL
    test_logical
    test_real
    test_float
    test_validated_float
    test_interval
    test_vector
    test_matrix
    test_series
    test_multi_index
    test_expansion
    test_polynomial
    test_differential
    test_procedure
    test_function
    test_expression
    test_taylor_model
    test_taylor_function
    test_function_sets
    test_affine_sets
    test_grid_set
    test_linear_programming
    test_nonlinear_programming
    test_integrator
    test_constraint_solver
    test_solvers
    test_differential_inclusion
    test_discrete_evolution
    test_continuous_evolution
    test_graphics
#    test_serialization
)

set(UNIT_TESTS_HYBRID
    test_hybrid_set
    test_hybrid_automaton
    test_hybrid_evolver
    test_hybrid_evolution
    test_hybrid_reachability_analyser
    test_crossings_issue
)

set(UNIT_TESTS ${UNIT_TESTS_KERNEL} ${UNIT_TESTS_HYBRID})

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})
endforeach()

add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})
add_custom_target(tests-kernel)
add_dependencies(tests-kernel ${UNIT_TESTS_KERNEL})
