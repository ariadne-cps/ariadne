add_custom_target(tests)
add_custom_target(tests-cpp)
add_dependencies(everything tests)

# Tests are grouped by architectural ownership rather than the former dynamics layer.
add_subdirectory(evolution)
add_subdirectory(analysis)
