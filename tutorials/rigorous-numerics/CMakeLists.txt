# Preferable to avoid warnings during execution
cmake_minimum_required (VERSION 3.1)

# Set the executable along with the required source files
add_executable(rigorous-numerics_tutorial rigorous-numerics_tutorial.cc)

# Instruct to link against the ariadne library
target_link_libraries(rigorous-numerics_tutorial ariadne)