# Preferable to avoid warnings during execution
cmake_minimum_required(VERSION 3.9)

project(RigorousNumericsTutorial VERSION 2.1)

# Set to search the root source directory for the FindAriadne.cmake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Set the C++ standard to C++20, as used by Ariadne
set(CMAKE_CXX_STANDARD 20)

# Find Ariadne and other dynamic library dependencies
find_package(MPFR REQUIRED)
find_package(Ariadne REQUIRED)

# Add the headers to the includes
include_directories(SYSTEM ${MPFR_INCLUDES} ${ARIADNE_INCLUDE_DIRS})

# Set the executable along with the required source files
add_executable(rigorous_numerics_tutorial rigorous_numerics_tutorial.cpp)

# Instruct to link against the ariadne library
target_link_libraries(rigorous_numerics_tutorial ${ARIADNE_LIBRARIES})