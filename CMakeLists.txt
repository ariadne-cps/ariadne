cmake_minimum_required (VERSION 2.6)
project (Ariadne)

# enable dashboard scripting
enable_testing()
include (CTest)

set (Ariadne_VERSION_MAJOR 0)
set (Ariadne_VERSION_MINOR 5)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include (CheckTypeSize)
include (CheckIncludeFile)
include (CheckIncludeFileCXX)
include (CheckLibraryExists)

check_type_size("int" SIZEOF_INT)

check_include_file (gmp.h HAVE_GMP_H)
check_include_file_cxx (gmpxx.h HAVE_GMPXX_H)
if(NOT HAVE_GMPXX_H)
    message(WARNING "Unable to find gmpxx.h.")
endif()

check_include_file (mpfr.h HAVE_MPFR_H)

check_include_file (cairo/cairo.h HAVE_CAIRO_H)
if(NOT HAVE_CAIRO_H)
    message(WARNING "Unable to find cairo.h.")
endif()

set (Boost_REQUIRED_VERSION 1.36.0)
find_package(Boost ${Boost_REQUIRED_VERSION})

IF(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
  set(Boost_USE_STATIC_LIBS OFF)
  set(Boost_USE_MULTITHREADED ON)
  set(Boost_USE_STATIC_RUNTIME OFF)
ELSEIF(NOT Boost_FOUND)
  message(FATAL_ERROR "Unable to find correct Boost version.")
ENDIF()

find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS serialization)
IF(Boost_SERIALIZATION_FOUND)
    set(HAVE_LIBBOOST_SERIALIZATION ON)
    set(ARIADNE_ENABLE_SERIALIZATION ON)
ELSEIF(NOT Boost_FOUND)
  message(FATAL_ERROR "Unable to find correct Boost Serialization version.")
ENDIF()

configure_file (${PROJECT_SOURCE_DIR}/config.h.cmake ${PROJECT_SOURCE_DIR}/config.h)

#Non-portable fix to avoid user header directories being searched for by system headers
add_definitions(-iquote${PROJECT_SOURCE_DIR})
#Should ideally be able to use line below
#include_directories(${PROJECT_SOURCE_DIR})


add_subdirectory(utility)
add_subdirectory(numeric)
add_subdirectory(algebra)
add_subdirectory(function)
add_subdirectory(solvers)
add_subdirectory(geometry)
add_subdirectory(dynamics)
add_subdirectory(expression)
add_subdirectory(output)
add_subdirectory(hybrid)

add_subdirectory(test)
add_subdirectory(wrap)
add_subdirectory(examples)
add_subdirectory(python)

add_library(ariadne SHARED
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-expression>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-dynamics>
    $<TARGET_OBJECTS:ariadne-output>
    $<TARGET_OBJECTS:ariadne-hybrid>
)

target_link_libraries(ariadne gmpxx gmp)
#FIXME: Make boost_system non-required
target_link_libraries(ariadne boost_system)
#FIXME: Make Cairo and Boost Serialization optional dependencies
target_link_libraries(ariadne cairo)
target_link_libraries (ariadne boost_serialization)

add_library(pyariadne $<TARGET_OBJECTS:pyariadne-src>)
set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne")
target_link_libraries (pyariadne ariadne boost_python)


add_library(pyinterval $<TARGET_OBJECTS:pyinterval-src>)
set_target_properties(pyinterval PROPERTIES PREFIX "" OUTPUT_NAME "interval")
target_link_libraries (pyinterval ariadne boost_python)

add_library(pydifferentiation $<TARGET_OBJECTS:pydifferentiation-src>)
set_target_properties(pydifferentiation PROPERTIES PREFIX "" OUTPUT_NAME "differentiation")
target_link_libraries (pydifferentiation ariadne boost_python)

