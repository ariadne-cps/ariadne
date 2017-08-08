set(EXAMPLES
    rigorous_numerics_tutorial
    hybrid_system_tutorial
    ballinabox
    bouncingball
    heating
    rectifier
)


foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(watertank)
add_subdirectory(twowatertanks)
add_subdirectory(experimental)
