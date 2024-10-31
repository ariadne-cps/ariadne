set(EXAMPLES_EXPERIMENTAL_DISCRETE
    heading_control
    discrete_probability_map
    continuous_probability_map
)

foreach(example ${EXAMPLES_EXPERIMENTAL_DISCRETE})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples ${EXAMPLES_EXPERIMENTAL_DISCRETE})

