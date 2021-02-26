set(TUTORIALS_DEMONSTRATION
    numeric_demonstration
    algebra_demonstration
    function_demonstration
    geometry_demonstration
    solver_demonstration
)

foreach(demonstration ${TUTORIALS_DEMONSTRATION})
    add_executable(${demonstration} ${demonstration}.cpp)
    target_link_libraries(${demonstration} ariadne)
endforeach()

add_dependencies(tutorials ${TUTORIALS_DEMONSTRATION})
