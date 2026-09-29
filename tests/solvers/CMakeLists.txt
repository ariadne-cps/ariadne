set(UNIT_TESTS
    test_bounder
    test_constraint_solver
    test_integrator    
    test_linear_programming
    test_nonlinear_programming
    test_solvers
    test_smt_solver
    test_smt_neural_benchmarks
    test_smt_boolean
    test_smt_theory
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel)
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "solvers;kernel")
endforeach()

target_compile_definitions(
    test_smt_neural_benchmarks
    PRIVATE ARIADNE_SMT_BARR3_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/smt_barr3_full64.bin")

add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})

add_executable(benchmark_smt_barr3_verification benchmark_smt_barr3_verification.cpp)
target_link_libraries(benchmark_smt_barr3_verification ariadne-kernel)
target_compile_definitions(
    benchmark_smt_barr3_verification
    PRIVATE ARIADNE_SMT_BARR3_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/smt_barr3_full64.bin")


add_executable(benchmark_constraint_propagation benchmark_constraint_propagation.cpp)
target_link_libraries(benchmark_constraint_propagation ariadne-kernel)
