# Preferable to avoid warnings during execution
cmake_minimum_required (VERSION 3.1)

# Set to search the root source directory for any .cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# Set the C++ standard to C++17, as used by Ariadne
set(CMAKE_CXX_STANDARD 17)

# Find Ariadne, adding the include directory if present, and failing if not
find_package(Ariadne)
if(ARIADNE_FOUND)
    include_directories(${ARIADNE_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Unable to find Ariadne.")
endif()

# Set the executable along with the required source files
add_executable(hybrid-evolution_tutorial hybrid-evolution_tutorial.cpp)

# Instruct to link against the ariadne library
target_link_libraries(hybrid-evolution_tutorial ${ARIADNE_LIBRARIES})