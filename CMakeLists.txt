cmake_minimum_required(VERSION 3.9)

project(Ariadne VERSION 1.9)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  message(WARNING "Intel compiler not supported.")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  message(WARNING "Microsoft compiler not supported.")
endif()

set(COMPILER_WARN all extra pedantic sign-conversion shadow cast-qual disabled-optimization 
                  init-self missing-include-dirs sign-promo switch-default undef redundant-decls)
set(COMPILER_DONT_WARN strict-aliasing unused-parameter)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_MACOSX_RPATH 1)
set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
foreach(WARN ${COMPILER_WARN})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W${WARN}")       
endforeach()
foreach(WARN ${COMPILER_DONT_WARN})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-${WARN}")       
endforeach()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-non-template-friend")
endif()
if(WERROR)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")
endif()

set(COVERAGE_COMPILER_FLAGS -O0 --coverage -fprofile-arcs -ftest-coverage)

enable_testing()
include(CTest)

include(cotire)

include(UseDoxygen OPTIONAL)

include(CheckTypeSize)
include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

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
    set(ARIADNE_ENABLE_STACK_TRACE ON)
endif()

find_package(GMP)
if(GMP_FOUND)
    include_directories(SYSTEM ${GMP_INCLUDES})
else()
    message(FATAL_ERROR "Unable to find GMP.")
endif()

find_package(GMPXX)
if(GMPXX_FOUND)
    set(HAVE_GMPXX_H ON)
    include_directories(SYSTEM ${GMPXX_INCLUDE_DIR})
else()
    message(WARNING "Unable to find GMPXX. Related functionality will be disabled.")
endif()    

find_package(MPFR)
if(MPFR_FOUND)
    include_directories(SYSTEM ${MPFR_INCLUDES})
else()
    message(FATAL_ERROR "Unable to find MPFR.")
endif()

find_package(Cairo)
if(CAIRO_FOUND)
    set(HAVE_CAIRO_H ON)
    include_directories(SYSTEM ${CAIRO_INCLUDE_DIRS})
else()
    message(WARNING "Unable to find Cairo. You will not be able to use graphics.")
endif()

find_package(GTK2)
if(GTK2_FOUND)
    set(HAVE_GTK2_H ON)
    include_directories(SYSTEM ${GTK2_INCLUDE_DIRS})
else()
    message(WARNING "Unable to find GTK2. You will not be able to use graphics windows.")
endif()


configure_file(${PROJECT_SOURCE_DIR}/cmake/config.hpp.cmake ${PROJECT_SOURCE_DIR}/source/config.hpp)

include_directories(${PROJECT_SOURCE_DIR}/source)

add_subdirectory(source)
add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(tutorials)
add_subdirectory(experimental EXCLUDE_FROM_ALL)

find_package(PythonLibs)
if(PYTHONLIBS_FOUND)
    message(STATUS "PythonLibs found: Python bindings will be built.")
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/pybind11/CMakeLists.txt)
        message(STATUS "pybind11 dependency not available.")
        find_package(Git)
        if(GIT_FOUND)
            message(STATUS "Downloading it...")
            execute_process(COMMAND git submodule update --init --recursive OUTPUT_QUIET ERROR_QUIET)
            message(STATUS "Downloaded pybind11 dependency successfully.")
        else()
            message(ERROR "Git not found, pybind11 dependency could not be downloaded.")
        endif()
    endif()
    add_subdirectory(external/pybind11)
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
)

add_library(ariadne SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-output>
    $<TARGET_OBJECTS:ariadne-symbolic>
    $<TARGET_OBJECTS:ariadne-dynamics>
    $<TARGET_OBJECTS:ariadne-hybrid>
)

if(COVERAGE)
    target_link_libraries(ariadne gcov)
    target_link_libraries(ariadne-kernel gcov)
    target_link_libraries(ariadne-core gcov)
endif()

target_link_libraries(ariadne dl)
target_link_libraries(ariadne-kernel dl)
target_link_libraries(ariadne-core dl)

target_link_libraries(ariadne ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
target_link_libraries(ariadne-kernel ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
target_link_libraries(ariadne-core ${GMP_LIBRARIES} ${MPFR_LIBRARIES})

if(GMPXX_FOUND)
    target_link_libraries(ariadne ${GMPXX_LIBRARIES})
    target_link_libraries(ariadne-kernel ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
    target_link_libraries(ariadne-core ${GMP_LIBRARIES} ${MPFR_LIBRARIES})    
endif()

if(CAIRO_FOUND)
    target_link_libraries(ariadne ${CAIRO_LIBRARIES})
    target_link_libraries(ariadne-kernel ${CAIRO_LIBRARIES})
endif()

if(GTK2_FOUND)
    target_link_libraries(ariadne ${GTK2_LIBRARIES})
    target_link_libraries(ariadne-kernel ${GTK2_LIBRARIES})
endif()

install(TARGETS ariadne DESTINATION lib)
install(DIRECTORY source/ DESTINATION include/ariadne FILES_MATCHING PATTERN "*.hpp")
