set(EXAMPLES_NOISY
    harmonic-oscillator
    clock
	vanderpol-noisy
	fitzhugh-nagumo
	laub-loomis
	vinograd
	sphere
	order7
	wiggins1873
)

foreach(example ${EXAMPLES_NOISY})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()
