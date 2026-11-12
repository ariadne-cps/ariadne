cmake_minimum_required(VERSION 3.15)
project(Ariadne VERSION 3.0)

# Bootstrap the shared CMake helpers before loading Reasoning.
if(NOT TARGET ariadne-build)
    find_package(Git REQUIRED)
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env "LC_ALL=C"
                "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}"
                submodule update --init -- ss/bui
        RESULT_VARIABLE _status)
    if(_status)
        message(FATAL_ERROR "Cannot initialize Build submodule")
    endif()
    include("${PROJECT_SOURCE_DIR}/ss/bui/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Build CHAINS
        "build" "reasoning|build")
endif()

setup_project()
add_subdirectory(ss)

# Reuse Reasoning's complete object set and public header dependency graph.
get_target_property(_objects reasoning ARIADNE_PROJECT_OBJECT_TARGETS)
setup_project_library(TARGET ariadne OBJECTS ${_objects} LINK_LIBRARIES reasoning)
# The shared helper creates the library only for standalone builds.
# Examples includes Ariadne as a subproject and needs the same concrete target.
if(NOT TARGET ariadne)
    add_library(ariadne ${LIBRARY_KIND})
    foreach(_object IN LISTS _objects)
        target_sources(ariadne PRIVATE $<TARGET_OBJECTS:${_object}>)
    endforeach()
    target_link_libraries(ariadne PUBLIC reasoning)
    install(TARGETS ariadne
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)
endif()
target_include_directories(ariadne PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne>)
ariadne_register_public_headers(TARGET ariadne
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne.hpp")
ariadne_install_dependency_bundle(TARGET ariadne DESTINATION include/ariadne)

# The aggregate Python extension is configured separately from the C++ library.
if(BUILD_PYTHON_BINDINGS AND TARGET pyariadne-reasoning)
    add_subdirectory(python)
endif()
