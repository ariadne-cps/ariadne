add_custom_target(tests)
add_custom_target(tests-cpp)

add_dependencies(everything tests)

add_subdirectory(function)
add_subdirectory(symbolic)
add_subdirectory(geometry)
add_subdirectory(solving)
add_subdirectory(io)
add_subdirectory(dynamics)
add_subdirectory(hybrid)