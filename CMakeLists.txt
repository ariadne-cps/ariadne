cmake_minimum_required(VERSION 3.15)
message(STATUS "CMake version ${CMAKE_VERSION}")

project(Ariadne VERSION 2.5.3)

include(CMakePackageConfigHelpers)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(WARNING "Intel compiler not supported.")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    message(WARNING "Microsoft compiler not supported.")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0)
        message(FATAL_ERROR "GCC version must be at least 12!")
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if ((NOT APPLE) AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16.0))
        message(FATAL_ERROR "Clang version must be at least 16!")
    endif()
else()
    message(WARNING "You are using an unidentified compiler! Compilation has only been tested with Clang and GCC.")
endif()

set(ANY_TARGET_WARN all extra pedantic sign-conversion cast-qual disabled-optimization
                  init-self missing-include-dirs sign-promo switch-default undef redundant-decls unused-parameter)
set(ANY_TARGET_DONT_WARN strict-aliasing)
set(LIBRARY_EXCLUSIVE_WARN shadow)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_MACOSX_RPATH 1)
if(APPLE)
    # Build-tree executables must load the just-built Ariadne dylibs rather
    # than an older installation from CMAKE_INSTALL_PREFIX (e.g. /usr/local).
    set(CMAKE_BUILD_RPATH "${PROJECT_BINARY_DIR}")
endif()
set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

foreach(WARN ${ANY_TARGET_WARN})
    add_compile_options(-W${WARN})
endforeach()
foreach(WARN ${ANY_TARGET_DONT_WARN})
    add_compile_options(-Wno-${WARN})
endforeach()
if(WARNINGS_AS_ERRORS)
    add_compile_options(-Werror)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wno-non-template-friend)
endif()
if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (NOT APPLE))
    add_compile_options(-frounding-math)
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-fdiagnostics-color=always)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-fcolor-diagnostics)
endif ()

option(COVERAGE "Enable coverage reporting" OFF)

enable_testing()
include(CTest)

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

list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/submodules/configuration/cmake")
include("${PROJECT_SOURCE_DIR}/submodules/configuration/cmake/ProjectOptions.cmake")

require_initialized_submodules()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        REFERENCES
            "${PROJECT_SOURCE_DIR}|submodules/configuration"
            "${PROJECT_SOURCE_DIR}/submodules/threading|submodules/configuration"
    )
endif()

if(COVERAGE)
    if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
        enable_project_coverage()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        include(GCCCodeCoverage)
    else()
        message(FATAL_ERROR
            "Coverage is supported only with GCC or macOS AppleClang/Clang.")
    endif()
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
add_custom_target(tests)

configure_file(${PROJECT_SOURCE_DIR}/cmake/config.hpp.cmake ${PROJECT_SOURCE_DIR}/source/config.hpp)

add_subdirectory(submodules)

include_directories(${PROJECT_SOURCE_DIR}/source)

add_subdirectory(source)
add_subdirectory(python)

add_subdirectory(tests EXCLUDE_FROM_ALL)
add_subdirectory(examples EXCLUDE_FROM_ALL)
add_subdirectory(tutorials EXCLUDE_FROM_ALL)
add_subdirectory(experimental EXCLUDE_FROM_ALL)

add_library(ariadne-core SHARED
    $<TARGET_OBJECTS:ariadne-foundations>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:THREADING_SRC>
    $<TARGET_OBJECTS:LOGGING_SRC>
    $<TARGET_OBJECTS:UTILITY_SRC>
    source/geometry/interval.cpp
)

add_library(ariadne-kernel SHARED
    $<TARGET_OBJECTS:ariadne-foundations>
    $<TARGET_OBJECTS:ariadne-numeric>
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
    $<TARGET_OBJECTS:ariadne-foundations>
    $<TARGET_OBJECTS:ariadne-numeric>
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
target_link_libraries(ariadne PUBLIC GMP::GMP MPFR::MPFR $<BUILD_INTERFACE:threading>)
target_link_libraries(ariadne PRIVATE dl ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ariadne-kernel PUBLIC $<BUILD_INTERFACE:threading> PRIVATE dl ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ariadne-core PUBLIC $<BUILD_INTERFACE:threading> PRIVATE dl ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

if(APPLE)
    # Keep the dylib identity relocatable.  In the build tree, CMake's
    # BUILD_RPATH above resolves @rpath to PROJECT_BINARY_DIR; installed
    # consumers can supply their own install/runtime rpath independently.
    set_target_properties(ariadne ariadne-kernel ariadne-core PROPERTIES
        MACOSX_RPATH ON
        INSTALL_NAME_DIR "@rpath"
        BUILD_WITH_INSTALL_NAME_DIR OFF
    )
endif()

if(COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        # Ariadne deliberately instruments only the library code on GCC.
        # Enabling coverage before the subdirectories would also instrument
        # tests and Python bindings, changing the gcov dataset.
        append_gcc_coverage_compiler_flags()
    endif()

    set(_ariadne_coverage_filter_args)
    if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
        list(APPEND _ariadne_coverage_filter_args
            EXCLUDE_REGEX ".*/(tests|python|submodules)/.*")
    endif()

    setup_project_coverage(
        TARGET ariadne
        TEST_TARGET tests
        SOURCES "${PROJECT_SOURCE_DIR}/source"
        ${_ariadne_coverage_filter_args}
    )
endif()

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
get_target_property(LOGGING_SOURCE_DIR logging SOURCE_DIR)
get_target_property(THREADING_SOURCE_DIR threading SOURCE_DIR)
get_target_property(UTILITY_SOURCE_DIR utility SOURCE_DIR)

install(DIRECTORY source/ DESTINATION include/ariadne FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY "${LOGGING_SOURCE_DIR}/include/logging" DESTINATION include/ariadne)
install(DIRECTORY "${THREADING_SOURCE_DIR}/include/threading" DESTINATION include/ariadne)
install(DIRECTORY "${UTILITY_SOURCE_DIR}/include/utility" DESTINATION include/ariadne)
