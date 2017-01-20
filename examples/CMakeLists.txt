set(EXAMPLES 
    tutorial
    ballinabox
    bouncingball
    watertank
    heating
    rectifier
    watertank-compositional
)


foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cc)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(experimental)