set(EXPERIMENTAL_EXAMPLES
	laser
	vanderpol-vectorfield
	vanderpol-nonoise
	harmonic-oscillator-nonoise
	lotka-volterra-nonoise
	higgins-selkov-nonoise
	laub-loomis-nonoise
)

foreach(example ${EXPERIMENTAL_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples-experimental ${EXPERIMENTAL_EXAMPLES})
