set(UNIT_TESTS
    test_bounder
    test_constraint_solver
    test_integrator    
    test_linear_programming
    test_nonlinear_programming
    test_solvers
    test_smt_solver
    test_smt_boolean
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel)
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "solvers;kernel")
endforeach()

add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})