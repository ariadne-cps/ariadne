set(EXPERIMENTAL_EXAMPLES
	circle
	harmonic-oscillator
	exponential
	laser
	bouncing-walls
	vanderpol
	vanderpol-vectorfield
	lotka-volterra
	higgins-selkov
	laub-loomis
)

foreach(example ${EXPERIMENTAL_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples-experimental ${EXPERIMENTAL_EXAMPLES})
