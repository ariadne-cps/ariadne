set(UNIT_TESTS_NUMERIC
    test_logical
    test_integer
    test_dyadic
    test_rational
    test_real
    test_number
    test_float
    test_validated_float
    test_rounding
)

set(UNIT_TESTS_ALGEBRA
    test_vector
    test_matrix
    test_series
    test_multi_index
    test_expansion
    test_differential
)

set(UNIT_TESTS_FUNCTION
    test_polynomial
    test_procedure
    test_function
    test_taylor_model
    test_taylor_function
        test_expression
)

set(UNIT_TESTS_GEOMETRY
    test_interval
    test_box
    test_function_sets
    test_affine_sets
    test_grid_set
)

set(UNIT_TESTS_SOLVERS
    test_linear_programming
    test_nonlinear_programming
    test_integrator
    test_constraint_solver
    test_solvers
)

set(UNIT_TESTS_DYNAMICS
    test_differential_inclusion
    test_discrete_evolution
    test_continuous_evolution
)

set(UNIT_TESTS_OUTPUT
    test_graphics
#    test_serialization
)

set(UNIT_TESTS_KERNEL ${UNIT_TESTS_NUMERIC} ${UNIT_TESTS_ALGEBRA} ${UNIT_TESTS_FUNCTION} ${UNIT_TESTS_GEOMETRY}
                      ${UNIT_TESTS_SOLVERS} ${UNIT_TESTS_DYNAMICS} ${UNIT_TESTS_OUTPUT})

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
    add_test(${TEST} ${TEST})
endforeach()

foreach(TEST ${UNIT_TESTS_NUMERIC} ${UNIT_TESTS_ALGEBRA})
    target_link_libraries(${TEST} ariadne-core)
endforeach()

foreach(TEST ${UNIT_TESTS_FUNCTION} ${UNIT_TESTS_GEOMETRY})
    target_link_libraries(${TEST} ariadne-kernel)
endforeach()

foreach(TEST ${UNIT_TESTS_SOLVERS} ${UNIT_TESTS_DYNAMICS} ${UNIT_TESTS_OUTPUT})
    target_link_libraries(${TEST} ariadne-kernel)
endforeach()

foreach(TEST ${UNIT_TESTS_HYBRID})
    target_link_libraries(${TEST} ariadne)
endforeach()

add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})
add_custom_target(tests-kernel)
add_dependencies(tests-kernel ${UNIT_TESTS_KERNEL})
