add_executable(
    benchmark_smt_barr3_verification
    benchmark_smt_barr3_verification.cpp)
target_link_libraries(
    benchmark_smt_barr3_verification
    ariadne-kernel)
target_compile_definitions(
    benchmark_smt_barr3_verification
    PRIVATE ARIADNE_SMT_BARR3_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/smt_barr3_full64.bin")

add_executable(
    benchmark_constraint_propagation
    benchmark_constraint_propagation.cpp)
target_link_libraries(
    benchmark_constraint_propagation
    ariadne-kernel)

add_dependencies(
    benchmarks
    benchmark_smt_barr3_verification
    benchmark_constraint_propagation)
