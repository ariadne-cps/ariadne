cmake_minimum_required(VERSION 3.15)
project(Ariadne VERSION 2.5.3)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

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
    require_same_dependency_commit(Configuration CHAINS
        "configuration"
        "evolution|configuration"
        "evolution|integration|configuration"
        "evolution|equation-solving|configuration"
        "evolution|system|configuration"
        "evolution|system|hybrid|configuration"
        "smt|configuration"
        "sat|configuration"
    )
    require_same_dependency_commit(Symbolic CHAINS
        "evolution|system|hybrid|labelled|expression-set|symbolic"
        "smt|symbolic"
    )
    require_same_dependency_commit(EquationSolving CHAINS
        "evolution|equation-solving"
        "smt|equation-solving"
    )
    require_same_dependency_commit(Function CHAINS
        "evolution|integration|function"
        "evolution|equation-solving|function"
        "evolution|system|hybrid|labelled|expression-set|symbolic|function"
    )
    require_same_dependency_commit(Cdcl CHAINS
        "smt|cdcl"
        "sat|cdcl"
    )
    require_same_dependency_commit(Bdd CHAINS
        "sat|bdd"
        "evolution|system|hybrid|labelled|drawing|paving|bdd"
    )
    require_same_dependency_commit(Threading CHAINS
        "smt|threading"
        "sat|bdd|threading"
    )
    require_same_dependency_commit(ConstraintSolving CHAINS
        "smt|constraint-solving"
        "evolution|system|hybrid|labelled|expression-set|set|constraint-solving"
    )
endif()

setup_project()

include_directories("${PROJECT_SOURCE_DIR}/source")
add_subdirectory(ss)
add_subdirectory(source)

# Reasoning is the only local object module; Evolution and its dependencies
# are supplied by the external Evolution submodule.
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
    $<TARGET_OBJECTS:ariadne-reasoning>
)

# Collect dependency objects exactly once, preserving the standalone
# repositories' aggregate interface and avoiding duplicate object linkage.
set(_ariadne_external_objects)
foreach(_dependency evolution smt sat)
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
    $<BUILD_INTERFACE:ariadne-reasoning>
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

# Local Reasoning bindings augment upstream Python registrations.
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

# Legacy examples/tutorials/experimental CMake references are excluded from
# configuration until their old solving/io headers have been migrated.
# Packaging must likewise be rebuilt around the new repository interfaces.
