cmake_minimum_required (VERSION 2.6)
project (Ariadne)

# enable dashboard scripting
enable_testing()
include (CTest)

set (Ariadne_VERSION_MAJOR 0)
set (Ariadne_VERSION_MINOR 5)

set(BUILD_SHARED_LIBS ON)

#configure_file ("${PROJECT_SOURCE_DIR}/config.h.in", "${PROJECT_BINARY_DIR}/config.h")

include_directories(${PROJECT_SOURCE_DIR})
# Non-portable fix to avoid user header directories being searched for by system headers
add_definitions(-iquote${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(wrap)
add_subdirectory(examples)
 

