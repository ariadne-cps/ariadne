cmake_minimum_required(VERSION 3.15)
message(STATUS "CMake version ${CMAKE_VERSION}")

project(Ariadne VERSION 2.5.3)

include(CMakePackageConfigHelpers)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

include(UseDoxygen OPTIONAL)

include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

check_include_file("dlfcn.h" HAVE_DLFCN_H)
check_include_file("execinfo.h" HAVE_EXECINFO_H)
check_include_file_cxx("cxxabi.h" HAVE_CXXABI_H)
if((${HAVE_DLFCN_H}) AND (${HAVE_EXECINFO_H}) AND (${HAVE_CXXABI_H}))
    set(ARIADNE_ENABLE_STACK_TRACE ON)
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/submodules/configuration/cmake/ProjectOptions.cmake")
    message(FATAL_ERROR
        "Git submodules are not initialized. "
        "Run 'git submodule update --init --recursive' before configuring the project.")
endif()

include("${PROJECT_SOURCE_DIR}/submodules/configuration/cmake/ProjectOptions.cmake")

require_initialized_submodules()

setup_project()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT APPLE)
    add_compile_options(-frounding-math)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        REFERENCES
            "${PROJECT_SOURCE_DIR}|submodules/configuration"
            "${PROJECT_SOURCE_DIR}/submodules/threading|submodules/configuration"
            "${PROJECT_SOURCE_DIR}/submodules/foundation|submodules/configuration"
            "${PROJECT_SOURCE_DIR}/submodules/numeric|submodules/configuration"
    )
    require_same_dependency_commit(Foundation
        REFERENCES
            "${PROJECT_SOURCE_DIR}|submodules/foundation"
            "${PROJECT_SOURCE_DIR}/submodules/numeric|submodules/foundation"
    )
    require_same_dependency_commit(Utility
        REFERENCES
            "${PROJECT_SOURCE_DIR}/submodules/threading/submodules/logging|submodules/utility"
            "${PROJECT_SOURCE_DIR}/submodules/foundation|submodules/utility"
    )
endif()


find_package(Threads)
if(NOT Threads_FOUND)
    message(FATAL_ERROR "Unable to find Threads library")
endif()

find_package(GMP)
if(GMP_FOUND)
    include_directories(SYSTEM ${GMP_INCLUDES})
else()
    message(FATAL_ERROR "Unable to find GMP")
endif()

find_package(MPFR)
if(MPFR_FOUND)
    include_directories(SYSTEM ${MPFR_INCLUDES})
else()
    message(FATAL_ERROR "Unable to find MPFR")
endif()

find_package(Cairo)
if(CAIRO_FOUND)
    set(HAVE_CAIRO_H ON)
    include_directories(SYSTEM ${CAIRO_INCLUDE_DIRS})
endif()

find_package(Gnuplot 5.2.8)
if(GNUPLOT_FOUND)
    set(HAVE_GNUPLOT_H ON)
endif()

if(NOT CAIRO_FOUND AND NOT GNUPLOT_FOUND)
    message(STATUS "No Cairo or Gnuplot backend support found: figures can't be written to file")
endif()

add_custom_target(everything)

configure_file(${PROJECT_SOURCE_DIR}/cmake/config.hpp.cmake ${PROJECT_SOURCE_DIR}/source/config.hpp)

add_subdirectory(submodules)

include_directories(${PROJECT_SOURCE_DIR}/source)

add_subdirectory(source)

add_subdirectory(examples EXCLUDE_FROM_ALL)
add_subdirectory(tutorials EXCLUDE_FROM_ALL)
add_subdirectory(experimental EXCLUDE_FROM_ALL)

add_library(ariadne-core SHARED
    $<TARGET_OBJECTS:FOUNDATION_SRC>
    $<TARGET_OBJECTS:NUMERIC_SRC>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:THREADING_SRC>
    $<TARGET_OBJECTS:LOGGING_SRC>
    $<TARGET_OBJECTS:UTILITY_SRC>
    source/geometry/interval.cpp
)

add_library(ariadne-kernel SHARED
    $<TARGET_OBJECTS:FOUNDATION_SRC>
    $<TARGET_OBJECTS:NUMERIC_SRC>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-io>
    $<TARGET_OBJECTS:ariadne-symbolic>
    $<TARGET_OBJECTS:THREADING_SRC>
    $<TARGET_OBJECTS:LOGGING_SRC>
    $<TARGET_OBJECTS:UTILITY_SRC>
)

add_library(ariadne SHARED
    $<TARGET_OBJECTS:FOUNDATION_SRC>
    $<TARGET_OBJECTS:NUMERIC_SRC>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-io>
    $<TARGET_OBJECTS:ariadne-symbolic>
    $<TARGET_OBJECTS:ariadne-dynamics>
    $<TARGET_OBJECTS:ariadne-hybrid>
    $<TARGET_OBJECTS:THREADING_SRC>
    $<TARGET_OBJECTS:LOGGING_SRC>
    $<TARGET_OBJECTS:UTILITY_SRC>
)

target_include_directories(ariadne PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
    $<INSTALL_INTERFACE:include/ariadne>
)
target_link_libraries(ariadne PUBLIC GMP::GMP MPFR::MPFR $<BUILD_INTERFACE:threading> $<BUILD_INTERFACE:foundation> $<BUILD_INTERFACE:numeric>)
target_link_libraries(ariadne PRIVATE dl ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ariadne-kernel PUBLIC $<BUILD_INTERFACE:threading> $<BUILD_INTERFACE:foundation> $<BUILD_INTERFACE:numeric> PRIVATE dl ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ariadne-core PUBLIC $<BUILD_INTERFACE:threading> $<BUILD_INTERFACE:foundation> $<BUILD_INTERFACE:numeric> PRIVATE dl ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

set(_ariadne_coverage_filter_args)
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
    list(APPEND _ariadne_coverage_filter_args
        EXCLUDE_REGEX ".*/(tests|python|submodules)/.*")
endif()

setup_standalone_project_tests(
    EXCLUDE_FROM_ALL
    TARGET ariadne
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES "${PROJECT_SOURCE_DIR}/source"
    ${_ariadne_coverage_filter_args}
)

add_subdirectory(python)

if(CAIRO_FOUND)
    target_link_libraries(ariadne PUBLIC Cairo::Cairo)
    target_link_libraries(ariadne-kernel PRIVATE ${CAIRO_LIBRARIES})
endif()

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/AriadneConfig.cmake.in
    ${PROJECT_BINARY_DIR}/AriadneConfig.cmake
    INSTALL_DESTINATION lib/cmake/Ariadne
)
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/AriadneConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(TARGETS ariadne EXPORT AriadneTargets DESTINATION lib)
install(EXPORT AriadneTargets
    FILE AriadneTargets.cmake
    NAMESPACE Ariadne::
    DESTINATION lib/cmake/Ariadne
)
install(FILES
    ${PROJECT_BINARY_DIR}/AriadneConfig.cmake
    ${PROJECT_BINARY_DIR}/AriadneConfigVersion.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindGMP.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindMPFR.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindCairo.cmake
    DESTINATION lib/cmake/Ariadne
)
get_target_property(FOUNDATION_SOURCE_DIR foundation SOURCE_DIR)
get_target_property(NUMERIC_SOURCE_DIR numeric SOURCE_DIR)
get_target_property(LOGGING_SOURCE_DIR logging SOURCE_DIR)
get_target_property(THREADING_SOURCE_DIR threading SOURCE_DIR)
get_target_property(UTILITY_SOURCE_DIR utility SOURCE_DIR)

install(DIRECTORY source/ DESTINATION include/ariadne FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY "${FOUNDATION_SOURCE_DIR}/include/foundation" DESTINATION include/ariadne)
install(DIRECTORY "${NUMERIC_SOURCE_DIR}/include/numeric" DESTINATION include/ariadne)
install(DIRECTORY "${LOGGING_SOURCE_DIR}/include/logging" DESTINATION include/ariadne)
install(DIRECTORY "${THREADING_SOURCE_DIR}/include/threading" DESTINATION include/ariadne)
install(DIRECTORY "${UTILITY_SOURCE_DIR}/include/utility" DESTINATION include/ariadne)
