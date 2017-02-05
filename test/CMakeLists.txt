set(UNIT_TESTS
    test_logical
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
    test_discrete_evolution
    test_continuous_evolution
    test_graphics
#    test_serialization
    test_hybrid_set
    test_hybrid_automaton
    test_hybrid_evolver
    test_hybrid_evolution
    test_hybrid_reachability_analyser
)

foreach(test ${UNIT_TESTS})
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} ariadne)
    add_test(${test} ${test})
endforeach()

