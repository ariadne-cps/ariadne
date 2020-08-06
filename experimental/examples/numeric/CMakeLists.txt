set(EXAMPLES_EXPERIMENTAL_NUMERIC
    real_calculus
    exact_real_computation
)

foreach(example ${EXAMPLES_EXPERIMENTAL_NUMERIC})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()
