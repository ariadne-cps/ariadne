set(EXAMPLES_EXPERIMENTAL_ARCH
    PRDE20
    CVDP20
	LALO20
    SUTR21
    LOVO21
    SPRE20
    arch_suite
)

foreach(example ${EXAMPLES_EXPERIMENTAL_ARCH})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()