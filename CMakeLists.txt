cmake_minimum_required (VERSION 3.1)
project (Ariadne)

set(Ariadne_VERSION_MAJOR 2)
set(Ariadne_VERSION_MINOR 0)

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  message(WARNING "Intel compiler not supported.")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  message(WARNING "Microsoft compiler not supported.")
endif()

enable_testing()
include(CTest)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_STANDARD 14)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_CXX_FLAGS "-O1")
endif()

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

check_include_file("dlfcn.h" HAVE_DLFCN_H)
check_include_file("execinfo.h" HAVE_EXECINFO_H)
check_include_file_cxx("cxxabi.h" HAVE_CXXABI_H)
if((${HAVE_DLFCN_H}) AND (${HAVE_EXECINFO_H}) AND (${HAVE_CXXABI_H}))
    set(ARIADNE_ENABLE_STACK_TRACE TRUE)
endif()

find_package(GMP)
if(GMP_FOUND)
    include_directories(${GMP_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Unable to find GMP.")
endif()

find_package(MPFR)
if(MPFR_FOUND)
    include_directories(${MPFR_INCLUDES})
else()
    message(FATAL_ERROR "Unable to find MPFR.")
endif()

find_package(Cairo)
if(CAIRO_FOUND)
    set(HAVE_CAIRO_H ON)
    include_directories( ${CAIRO_INCLUDE_DIRS})
else()
    message(WARNING "Unable to find Cairo. You will not be able to use graphics.")
endif()

find_package(GTK2)
if(GTK2_FOUND)
    set(HAVE_GTK2_H ON)
    include_directories( ${GTK2_INCLUDE_DIRS})
    link_directories( ${GTK2_LIBRARY_DIRS})
else()
    message(WARNING "Unable to find GTK2. You will not be able to use graphics windows.")
endif()

set(Boost_REQUIRED_VERSION 1.36.0)
find_package(Boost ${Boost_REQUIRED_VERSION})

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(Boost_USE_STATIC_LIBS OFF)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_RUNTIME OFF)
elseif(NOT Boost_FOUND)
    message(WARNING "Unable to find correct Boost version (>=${Boost_REQUIRED_VERSION})")
endif()

find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS system serialization)
if(Boost_SERIALIZATION_FOUND)
    set(HAVE_LIBBOOST_SERIALIZATION_HPP ON)
    set(ARIADNE_ENABLE_SERIALIZATION ON)
    set(BOOST_SERIALIZATION_LIBRARIES "${Boost_LIBRARIES}")
elseif(Boost_FOUND)
    message(WARNING "Unable to find Boost Serialization component.")
endif()

find_package(PythonLibs 2.4)
find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS python)
find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS python27)
if(PYTHONLIBS_FOUND AND (Boost_PYTHON_FOUND OR Boost_PYTHON27_FOUND))
    set(USE_PYTHON ON)
    set(BOOST_PYTHON_LIBRARIES "${Boost_LIBRARIES}") 
    message("Python libraries will be built")
elseif(NOT PYTHONLIBS_FOUND)
    message(WARNING "PythonLibs not found: cannot use Python bindings.")
elseif(NOT(Boost_PYTHON_FOUND OR Boost_PYTHON27_FOUND))
    message(WARNING "Unable to find any Boost Python component: cannot use Python bindings.")
endif()

configure_file (${PROJECT_SOURCE_DIR}/config.hpp.cmake ${PROJECT_SOURCE_DIR}/source/config.hpp)

include_directories(${PROJECT_SOURCE_DIR}/source)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  add_definitions(-Wno-non-template-friend)
endif()

add_subdirectory(source)
add_subdirectory(test)
add_subdirectory(examples)
add_subdirectory(tutorials)

if(USE_PYTHON)
    add_subdirectory(wrap)
    add_subdirectory(python)
endif()

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
    $<TARGET_OBJECTS:ariadne-symbolic>
    $<TARGET_OBJECTS:ariadne-dynamics>
)

add_library(ariadne SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-symbolic>
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

target_link_libraries(ariadne ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
target_link_libraries(ariadne-kernel ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
target_link_libraries(ariadne-core ${GMP_LIBRARIES} ${MPFR_LIBRARIES})

if(Boost_SERIALIZATION_FOUND)
    target_link_libraries(ariadne ${BOOST_SERIALIZATION_LIBRARIES})
endif()

if(GTK2_FOUND)
    target_link_libraries(ariadne ${GTK2_LIBRARIES})
endif()

if(CAIRO_FOUND)
    target_link_libraries(ariadne ${CAIRO_LIBRARIES})
    target_link_libraries(ariadne-kernel ${CAIRO_LIBRARIES})
endif()

install (TARGETS ariadne DESTINATION lib)
install (DIRECTORY source/ DESTINATION include/ariadne
            FILES_MATCHING PATTERN "*.hpp"
            PATTERN "source/prototyping" EXCLUDE
            PATTERN "source/experimental" EXCLUDE)

