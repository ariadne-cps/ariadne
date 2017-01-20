cmake_minimum_required (VERSION 2.8.12)
project (Ariadne)

if(POLICY CMP0045)
  cmake_policy(SET CMP0045 OLD)
endif()
if(POLICY CMP0058)
  cmake_policy(SET CMP0058 NEW)
endif()

# enable dashboard scripting
enable_testing()
include(CTest)

set(Ariadne_VERSION_MAJOR 0)
set(Ariadne_VERSION_MINOR 5)

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  message(WARNING "Intel compiler not supported.")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  message(WARNING "Microsoft compiler not supported.")
endif()

set(BUILD_SHARED_LIBS ON)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O1")

# add support for precompiled headers and single compilation unit
include(cotire)

include (CheckTypeSize)
include (CheckIncludeFile)
include (CheckIncludeFileCXX)
include (CheckLibraryExists)

check_type_size("int" SIZEOF_INT)

check_include_file(gmp.h HAVE_GMP_H)
check_include_file_cxx(gmpxx.h HAVE_GMPXX_H)
if(NOT HAVE_GMPXX_H)
    message(WARNING "Unable to find gmpxx.h.")
endif()

check_include_file(mpfr.h HAVE_MPFR_H)

find_package(GTK2)

check_include_file(cairo/cairo.h HAVE_CAIRO_H)
if(NOT HAVE_CAIRO_H)
    message(WARNING "Unable to find cairo.h.")
endif()

set(Boost_REQUIRED_VERSION 1.36.0)
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
#add_definitions(-iquote${PROJECT_SOURCE_DIR})
#Should ideally be able to use line below
include_directories(${PROJECT_SOURCE_DIR} ${CAIRO_INCLUDE_DIRS})

IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  add_definitions(-Wno-non-template-friend)
ENDIF()

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
add_subdirectory(python)
add_subdirectory(examples)

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
#FIXME: Make mpfr non-required
target_link_libraries(ariadne mpfr)
#FIXME: Make boost_system non-required
target_link_libraries(ariadne boost_system)
#FIXME: Make Cairo and Boost Serialization optional dependencies
target_link_libraries(ariadne cairo)
target_link_libraries (ariadne boost_serialization)

#add_library(pyariadne $<TARGET_OBJECTS:pyariadne-src>)
#set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne")
#target_link_libraries (pyariadne ariadne boost_python)

#add_library(pynumeric $<TARGET_OBJECTS:pynumeric-src>)
#set_target_properties(pynumeric PROPERTIES PREFIX "" OUTPUT_NAME "numeric")
#target_link_libraries (pynumeric ariadne boost_python)

#add_library(pyalgebra $<TARGET_OBJECTS:pyalgebra-src>)
#set_target_properties(pyalgebra PROPERTIES PREFIX "" OUTPUT_NAME "algebra")
#target_link_libraries (pyalgebra ariadne boost_python)
