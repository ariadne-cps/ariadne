set(EXAMPLES_EXPERIMENTAL_CONTINUOUS
    brusselator
    laub-loomis
    inner
)

foreach(example ${EXAMPLES_EXPERIMENTAL_CONTINUOUS})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(noisy)
add_subdirectory(constrained)