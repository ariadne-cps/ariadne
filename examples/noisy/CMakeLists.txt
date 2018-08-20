set(NOISY_EXAMPLES
	higgins-selkov
)

foreach(example ${NOISY_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples-noisy ${NOISY_EXAMPLES})
