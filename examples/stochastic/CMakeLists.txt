set(EXAMPLES_STOCHASTIC
    gaussian
)

foreach(example ${EXAMPLES_STOCHASTIC})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples ${EXAMPLES_STOCHASTIC})

