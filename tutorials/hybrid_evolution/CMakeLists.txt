# Preferable to avoid warnings during execution
cmake_minimum_required(VERSION 3.9)

# Set to search the root source directory for the FindAriadne.cmake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# Set to look for includes in usr/local directory, required for GMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include")

# Set the C++ standard to C++17, as used by Ariadne
set(CMAKE_CXX_STANDARD 17)

# Find Ariadne as a required dependency
find_package(Ariadne REQUIRED)

# Add the Ariadne headers to the includes
include_directories(SYSTEM ${ARIADNE_INCLUDE_DIRS})

# Set the executable along with the required source files
add_executable(hybrid_evolution_tutorial hybrid_evolution_tutorial.cpp)

# Instruct to link against the ariadne library
target_link_libraries(hybrid_evolution_tutorial ${ARIADNE_LIBRARIES})