set(UNIT_TESTS
    test_constraint_solver
    test_integrator    
    test_linear_programming
    test_nonlinear_programming
    test_solvers
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    add_test(${TEST} ${TEST})
    target_link_libraries(${TEST} ariadne-kernel)
endforeach()

add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})