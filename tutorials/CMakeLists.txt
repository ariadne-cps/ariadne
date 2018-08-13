set(TUTORIALS
    rigorous_numerics_tutorial
    hybrid_system_tutorial
)


foreach(tutorial ${TUTORIALS})
    add_executable(${tutorial} ${tutorial}.cpp)
    target_link_libraries(${tutorial} ariadne)
endforeach()

add_custom_target(tutorials)
add_dependencies(tutorials ${TUTORIALS})
