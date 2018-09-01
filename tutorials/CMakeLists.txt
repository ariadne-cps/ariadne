set(TUTORIALS
    rigorous-numerics
    hybrid-evolution
)

add_custom_target(tutorials)

foreach(tutorial ${TUTORIALS})
    add_executable(${tutorial}_tutorial ${tutorial}/${tutorial}_tutorial.cpp)
    target_link_libraries(${tutorial}_tutorial ariadne)
    add_dependencies(tutorials ${tutorial}_tutorial)
endforeach()

