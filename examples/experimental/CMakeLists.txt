set(EXPERIMENTAL_EXAMPLES
	circle
	exponential
	laser
	bouncing-walls
	vanderpol
)

foreach(example ${EXPERIMENTAL_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()
