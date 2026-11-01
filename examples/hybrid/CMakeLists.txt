set(EXAMPLES_HYBRID
    bouncingball
    heating
    rectifier
    laser
)

foreach(example ${EXAMPLES_HYBRID})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples ${EXAMPLES_HYBRID})

add_subdirectory(watertank)
add_subdirectory(twowatertanks)
add_subdirectory(threewatertanks)
