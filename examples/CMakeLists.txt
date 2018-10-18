set(EXAMPLES
    ballinabox
    bouncingball
    heating
    rectifier
)

foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_custom_target(examples)
add_dependencies(examples ${EXAMPLES})

add_subdirectory(watertank)
add_subdirectory(twowatertanks)
add_subdirectory(threewatertanks)
add_subdirectory(noisy)
add_subdirectory(pdes)
