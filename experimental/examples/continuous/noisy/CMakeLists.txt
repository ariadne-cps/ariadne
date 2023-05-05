set(EXAMPLES_EXPERIMENTAL_NOISY
    harmonic-oscillator
    clock
	vanderpol-noisy
	fitzhugh-nagumo
	vinograd
	sphere
	order7
	wiggins1873
	ross-mcdonald
	michelson
)

foreach(example ${EXAMPLES_EXPERIMENTAL_NOISY})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()
