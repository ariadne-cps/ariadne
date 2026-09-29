set(BENCHMARKS
    benchmark_taylor_model
)

foreach(BENCHMARK ${BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne)
endforeach()

add_dependencies(benchmarks ${BENCHMARKS})

add_subdirectory(solvers)
