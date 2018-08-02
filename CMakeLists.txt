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

include(UseDoxygen OPTIONAL)

include (CheckTypeSize)
include (CheckIncludeFile)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckLibraryExists)

check_type_size("char" SIZEOF_CHAR)
check_type_size("short int" SIZEOF_SHORT_INT)
check_type_size("int" SIZEOF_INT)
check_type_size("long int" SIZEOF_LONG_INT)
check_type_size("long long int" SIZEOF_LONG_LONG_INT)
check_type_size("float" SIZEOF_FLOAT)
check_type_size("double" SIZEOF_DOUBLE)
check_type_size("long double" SIZEOF_LONG_DOUBLE)

check_include_file(gmp.h HAVE_GMP_H)
check_include_file(mpfr.h HAVE_MPFR_H)
check_include_file_cxx(gmpxx.h HAVE_GMPXX_H)

check_include_file("dlfcn.h" HAVE_DLFCN_H)
check_include_file("execinfo.h" HAVE_EXECINFO_H)
check_include_file_cxx("cxxabi.h" HAVE_CXXABI_H)
if((${HAVE_DLFCN_H}) AND (${HAVE_EXECINFO_H}) AND (${HAVE_CXXABI_H}))
    set(ARIADNE_ENABLE_STACK_TRACE TRUE)
endif()

if(NOT HAVE_GMP_H)
    message(FATAL_ERROR "Unable to find required header gmp.h")
endif()

if(NOT HAVE_MPFR_H)
    message(FATAL_ERROR "Unable to find required header mpfr.h")
endif()

find_package(Cairo)
if(CAIRO_FOUND)
    message("Found drawing package Cairo")
    #message("  CAIRO_INCLUDE_DIRS: ${CAIRO_INCLUDE_DIRS}")
    #message("  CAIRO_LIBRARY_DIRS: ${CAIRO_LIBRARY_DIRS}")
    #message("  CAIRO_LIBRARIES: ${CAIRO_LIBRARIES}")
    set(HAVE_CAIRO_H ON)
    include_directories( ${CAIRO_INCLUDE_DIRS})
    link_directories( ${CAIRO_LIBRARY_DIRS})
else()
    message(WARNING "Unable to find Cairo. You will not be able to use graphics.")
endif()

find_package(GTK2)
if(GTK2_FOUND)
    message("Found GTK2")
    #message("  GTK2_INCLUDE_DIRS: ${GTK2_INCLUDE_DIRS}")
    #message("  GTK2_LIBRARY_DIRS: ${GTK2_LIBRARY_DIRS}")
    #message("  GTK2_LIBRARIES: ${GTK2_LIBRARIES}")
    set(HAVE_GTK2_H ON)
    include_directories( ${GTK2_INCLUDE_DIRS})
    link_directories( ${GTK2_LIBRARY_DIRS})
else()
    message(WARNING "Unable to find GTK2. You will not be able to use graphics windows.")
endif()

set(Boost_REQUIRED_VERSION 1.36.0)
find_package(Boost ${Boost_REQUIRED_VERSION})

if(Boost_FOUND)
    #message("Found Boost")
    #message("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    #message("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
    #message("  Boost_LIBRARIES: ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(Boost_USE_STATIC_LIBS OFF)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_RUNTIME OFF)
elseif(NOT Boost_FOUND)
    message(FATAL_ERROR "Unable to find correct Boost version (>=${Boost_REQUIRED_VERSION})")
endif()

find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS serialization)
if(Boost_SERIALIZATION_FOUND)
    set(HAVE_LIBBOOST_SERIALIZATION_HPP ON)
    set(ARIADNE_ENABLE_SERIALIZATION ON)
elseif(NOT Boost_FOUND)
    message(WARNING "Unable to find Boost Serialization component.")
endif()

configure_file (${PROJECT_SOURCE_DIR}/config.h.cmake ${PROJECT_SOURCE_DIR}/source/config.h)

#Non-portable fix to avoid user header directories being searched for by system headers
#add_definitions(-iquote${PROJECT_SOURCE_DIR})
#Should ideally be able to use line below
include_directories(${PROJECT_SOURCE_DIR}/source)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  add_definitions(-Wno-non-template-friend)
endif()

add_subdirectory(source)
add_subdirectory(test)
add_subdirectory(wrap)
add_subdirectory(python)
#add_subdirectory(examples)
#add_subdirectory(tutorials)

add_library(ariadne-core SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-algebra>
    source/geometry/interval.cpp
)

add_library(ariadne-kernel SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-output>
    $<TARGET_OBJECTS:ariadne-expression>
    $<TARGET_OBJECTS:ariadne-dynamics>
)

add_library(ariadne SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-expression>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-output>
    $<TARGET_OBJECTS:ariadne-dynamics>
    $<TARGET_OBJECTS:ariadne-hybrid>
)

target_link_libraries(ariadne dl)
target_link_libraries(ariadne-kernel dl)
target_link_libraries(ariadne-core dl)

target_link_libraries(ariadne gmp mpfr)
target_link_libraries(ariadne-kernel gmp mpfr)
target_link_libraries(ariadne-core gmp mpfr)

if(Boost_SERIALIZATION_FOUND)
    target_link_libraries(ariadne boost_system boost_serialization)
endif()

if(GTK2_FOUND)
    target_link_libraries(ariadne ${GTK2_LIBRARIES})
endif()

if(CAIRO_FOUND)
    target_link_libraries(ariadne cairo)
    target_link_libraries(ariadne-kernel cairo)
endif()

find_package(PythonLibs 2.4)

add_library(pyariadne $<TARGET_OBJECTS:pyariadne-src>)
set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne")
target_link_libraries (pyariadne ariadne ${PYTHON_LIBRARIES} boost_python)

add_library(pyariadne-numeric SHARED $<TARGET_OBJECTS:pyariadne-numeric-src> $<TARGET_OBJECTS:ariadne-utility> $<TARGET_OBJECTS:ariadne-numeric>)
set_target_properties(pyariadne-numeric PROPERTIES PREFIX "" OUTPUT_NAME "numeric")
target_link_libraries (pyariadne-numeric ${PYTHON_LIBRARIES} boost_python mpfr gmp)

#add_library(pyariadne-algebra $<TARGET_OBJECTS:pyariadne-algebra-src>)
#set_target_properties(pyariadne-algebra PROPERTIES PREFIX "" OUTPUT_NAME "algebra")
#target_link_libraries (pyariadne-algebra ariadne boost_python)

install (TARGETS ariadne DESTINATION lib)
install (FILES source/config.h DESTINATION include/ariadne)
install (DIRECTORY source/ DESTINATION include/ariadne/
            FILES_MATCHING PATTERN "*.hpp"
            PATTERN "source/prototyping" EXCLUDE
            PATTERN "source/experimental" EXCLUDE)

