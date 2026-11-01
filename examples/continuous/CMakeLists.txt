set(EXAMPLES_CONTINUOUS
    attractor
    vanderpol
    lorenz
    laub-loomis
)

foreach(example ${EXAMPLES_CONTINUOUS})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples ${EXAMPLES_CONTINUOUS})

add_subdirectory(noisy)

