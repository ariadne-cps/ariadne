set(EXAMPLES_EXPERIMENTAL
	laser
	vanderpol-vectorfield
	vanderpol-nonoise
	harmonic-oscillator-nonoise
	lotka-volterra-nonoise
	higgins-selkov-nonoise
	laub-loomis-nonoise
)

foreach(example ${EXAMPLES_EXPERIMENTAL})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_custom_target(examples-experimental)
add_dependencies(examples-experimental ${EXAMPLES_EXPERIMENTAL})
add_dependencies(examples ${EXAMPLES_EXPERIMENTAL})
