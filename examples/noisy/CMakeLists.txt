set(EXAMPLES_NOISY
	higgins-selkov
	chemical-reactor
	lotka-volterra
	jet-engine
	pi-controller
	jerk21
	lorenz-attractor
	rossler-attractor
	jerk16
	dc-dc
	harmonic-oscillator
	
	noisy-benchmark
	
	
)

foreach(example ${EXAMPLES_NOISY})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_custom_target(examples-noisy)
add_dependencies(examples-noisy ${EXAMPLES_NOISY})
