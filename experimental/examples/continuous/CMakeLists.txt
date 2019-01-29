set(EXAMPLES_EXPERIMENTAL_CONTINUOUS
    laub-loomis-nonoise
	lotka-volterra-nonoise
	higgins-selkov-nonoise
    quadrotor
)

foreach(example ${EXAMPLES_EXPERIMENTAL_CONTINUOUS})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()