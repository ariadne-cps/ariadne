set(EXAMPLES_EXPERIMENTAL_ARCH
    ROBE25
    CVDP20
	LALO20
    SUTR21
    LOVO25
    SPRE20
    arch_suite
)

foreach(example ${EXAMPLES_EXPERIMENTAL_ARCH})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()