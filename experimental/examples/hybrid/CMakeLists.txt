set(EXAMPLES_EXPERIMENTAL_HYBRID
	laser
    space_rendezvous
    power_converters
)

foreach(example ${EXAMPLES_EXPERIMENTAL_HYBRID})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(twowatertanks)
add_subdirectory(threewatertanks)