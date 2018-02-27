set(EXPERIMENTAL_EXAMPLES
	circle
	exponential
	laser
	bouncing-walls
	vanderpol
	vanderpol-vectorfield
	harmonic_inclusion
)

foreach(example ${EXPERIMENTAL_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples-experimental ${EXPERIMENTAL_EXAMPLES})
