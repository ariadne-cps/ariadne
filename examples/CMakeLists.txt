set(EXAMPLES 
    tutorial
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
add_subdirectory(experimental)