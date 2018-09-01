# Preferable to avoid warnings during execution
cmake_minimum_required (VERSION 3.1)

# Set the executable along with the required source files
add_executable(hybrid-evolution_tutorial hybrid-evolution_tutorial.cc)

# Instruct to link against the ariadne library
target_link_libraries(hybrid-evolution_tutorial ariadne)