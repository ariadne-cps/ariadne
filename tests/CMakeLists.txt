add_custom_target(tests)
add_custom_target(tests-cpp)
add_dependencies(everything tests)

# Evolution tests are owned by the external Evolution repository.
add_subdirectory(reasoning)
