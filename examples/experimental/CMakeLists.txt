set(EXPERIMENTAL_EXAMPLES
	circle
	exponential
	laser
)

foreach(example ${EXPERIMENTAL_EXAMPLES})
    add_executable(${example} ${example}.cc)
    target_link_libraries(${example} ariadne)
endforeach()
