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
    set_property(TEST ${TEST} PROPERTY LABELS "solving;kernel")
endforeach()

target_include_directories(
    test_smt_neural_benchmarks
    PRIVATE ${PROJECT_SOURCE_DIR}/benchmarks/solvers)
target_compile_definitions(
    test_smt_neural_benchmarks
    PRIVATE ARIADNE_SMT_BARR3_DATA_PATH="${PROJECT_SOURCE_DIR}/benchmarks/solving/data/smt_barr3_full64.bin")

add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
