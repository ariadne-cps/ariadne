set(EXAMPLES_EXPERIMENTAL_CONTINUOUS
    brusselator
    laub-loomis
    inner
    constrained_vanderpol
)

foreach(example ${EXAMPLES_EXPERIMENTAL_CONTINUOUS})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(noisy)