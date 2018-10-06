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