set(EXAMPLES_EXPERIMENTAL_ARCH
    CVDP20
	LALO20
    LOVO21
    PRDE20
    SPRE20
    arch_suite
)

foreach(example ${EXAMPLES_EXPERIMENTAL_ARCH})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()