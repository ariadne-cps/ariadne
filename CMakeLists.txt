cmake_minimum_required(VERSION 3.15)
project(Ariadne VERSION 2.5.3)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
# The old aggregate Python bindings are not compatible with this split yet.
option(BUILD_PYTHON_BINDINGS "Build Python bindings after layer migration" OFF)

find_package(Git REQUIRED)
if(NOT TARGET ariadne-configuration)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}"
                submodule update --init -- ss/cfg
        RESULT_VARIABLE _configuration_result
    )
    if(_configuration_result)
        message(FATAL_ERROR "Cannot initialize Configuration submodule")
    endif()
    include("${PROJECT_SOURCE_DIR}/ss/cfg/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Compare immutable gitlink commits, not branch names.
    require_same_dependency_commit(Configuration CHAINS
        "configuration"
        "integration|configuration"
        "equation-solving|configuration"
        "smt|configuration"
        "sat|configuration"
        "labelling|configuration"
        "symbolic|configuration"
    )
    require_same_dependency_commit(Symbolic CHAINS
        "symbolic"
        "smt|symbolic"
        "labelling|expression-set|symbolic"
    )
    require_same_dependency_commit(EquationSolving CHAINS
        "equation-solving"
        "smt|equation-solving"
    )
    require_same_dependency_commit(Function CHAINS
        "integration|function"
        "equation-solving|function"
        "symbolic|function"
    )
    require_same_dependency_commit(Cdcl CHAINS
        "smt|cdcl"
        "sat|cdcl"
    )
    require_same_dependency_commit(Bdd CHAINS
        "sat|bdd"
        "labelling|drawing|paving|bdd"
    )
    require_same_dependency_commit(Threading CHAINS
        "smt|threading"
        "sat|bdd|threading"
    )
    require_same_dependency_commit(ConstraintSolving CHAINS
        "smt|constraint-solving"
        "labelling|expression-set|set|constraint-solving"
    )
endif()

setup_project()

# Legacy dynamics/hybrid code still includes "config.hpp".
# Generate it in the build tree, not the source checkout.
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/generated")
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/config.hpp.cmake"
    "${PROJECT_BINARY_DIR}/generated/config.hpp"
)
include_directories(BEFORE "${PROJECT_BINARY_DIR}/generated")
include_directories("${PROJECT_SOURCE_DIR}/source")
add_subdirectory(ss)
add_subdirectory(source)

# System definitions and the two remaining legacy object modules are Ariadne-owned.
# All other code is
# provided by upstream submodules and their transitively exposed dependencies.
# These transitional modules do not annotate exported symbols. On Windows,
# a DLL with no exports does not produce the import library needed by the
# test executables; use a static archive instead of a non-exporting DLL.
# Preserve the shared library on Unix platforms.
if(WIN32)
    set(_ariadne_library_type STATIC)
else()
    set(_ariadne_library_type SHARED)
endif()
add_library(ariadne ${_ariadne_library_type}
    $<TARGET_OBJECTS:ariadne-system>
    $<TARGET_OBJECTS:ariadne-dynamics>
    $<TARGET_OBJECTS:ariadne-hybrid>
)

# Collect dependency objects exactly once, preserving the standalone
# repositories' aggregate interface and avoiding duplicate object linkage.
set(_ariadne_external_objects)
foreach(_dependency integration equation_solving smt sat labelling symbolic)
    get_target_property(_objects ${_dependency} ARIADNE_PROJECT_OBJECT_TARGETS)
    if(_objects)
        list(APPEND _ariadne_external_objects ${_objects})
    endif()
endforeach()
list(REMOVE_DUPLICATES _ariadne_external_objects)
foreach(_object IN LISTS _ariadne_external_objects)
    target_sources(ariadne PRIVATE $<TARGET_OBJECTS:${_object}>)
endforeach()

target_include_directories(ariadne PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
    $<INSTALL_INTERFACE:include/ariadne>
)
target_link_libraries(ariadne PUBLIC
    $<BUILD_INTERFACE:ariadne-analysis>
)
if(UNIX AND NOT APPLE)
    target_link_libraries(ariadne PRIVATE dl)
endif()
find_package(Threads REQUIRED)
target_link_libraries(ariadne PRIVATE Threads::Threads)

setup_standalone_project_tests(
    EXCLUDE_FROM_ALL
    TARGET ariadne
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES "${PROJECT_SOURCE_DIR}/source"
)

# External dependencies own their Python bindings. The legacy Ariadne Python
# aggregate cannot be re-enabled before its layer-specific sources are split.
if(BUILD_PYTHON_BINDINGS)
    message(FATAL_ERROR
        "The legacy pyariadne bindings require migration to system/evolution/analysis.")
endif()

# Legacy examples/tutorials/experimental CMake references are excluded from
# configuration until their old solving/io headers have been migrated.
# Packaging must likewise be rebuilt around the new repository interfaces.
