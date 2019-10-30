add_custom_target(examples)

add_subdirectory(discrete)
add_subdirectory(continuous)
add_subdirectory(hybrid)
add_subdirectory(pdes)

add_custom_target(examples-external)
add_subdirectory(external)
