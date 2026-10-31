add_custom_target(tests)
add_custom_target(tests-cpp)
add_dependencies(everything tests)

# Tests are grouped by architectural ownership rather than the former dynamics layer.
# Only three historical system-model regressions remain locally until their precise
# scenarios are covered by the standalone system repository.
add_subdirectory(evolution)
add_subdirectory(analysis)
add_subdirectory(system_legacy)
