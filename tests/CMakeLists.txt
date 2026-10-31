add_custom_target(tests)
add_custom_target(tests-cpp)
add_dependencies(everything tests)

# Retired local Function/Geometry/Solving/Symbolic/IO/Hybrid/System test modules are
# no longer registered. Their tests belong to the corresponding repositories.
add_subdirectory(dynamics)
add_subdirectory(evolution)
add_subdirectory(analysis)
