set(EXAMPLES_EXPERIMENTAL_ARCH
    ROBE21
    CVDP23
	LALO20
    LOVO21
    SPRE22
    arch_suite
)

foreach(example ${EXAMPLES_EXPERIMENTAL_ARCH})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()