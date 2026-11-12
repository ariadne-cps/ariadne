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
target_include_directories(ariadne PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne>)
ariadne_register_public_headers(TARGET ariadne
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne.hpp")
ariadne_install_dependency_bundle(TARGET ariadne DESTINATION include/ariadne)

# Compose the same registered Python bindings under the pyariadne entry point.
if(TARGET pyariadne-reasoning)
    get_property(_modules GLOBAL PROPERTY ARIADNE_PYTHON_MODULES)
    foreach(_module IN LISTS _modules)
        string(APPEND _declarations "void ${_module}_module(pybind11::module&);\n")
        string(APPEND _initializers "    ${_module}_module(module);\n")
    endforeach()
    file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pyariadne_entry.cpp"
        CONTENT "#include \"pybind11.hpp\"\n${_declarations}\nPYBIND11_MODULE(pyariadne, module) {\n${_initializers}}\n")
    pybind11_add_module(pyariadne "${CMAKE_CURRENT_BINARY_DIR}/pyariadne_entry.cpp")
    target_link_libraries(pyariadne PRIVATE ariadne pyariadne-reasoning pyariadne-python-common)
    foreach(_object IN LISTS _objects)
        target_sources(pyariadne PRIVATE $<TARGET_OBJECTS:${_object}>)
    endforeach()
    if(WIN32)
        set(_python_install_dir "Lib/site-packages")
    else()
        set(_python_install_dir "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
    endif()
    install(TARGETS pyariadne DESTINATION "${_python_install_dir}")
    add_dependencies(bindings pyariadne)
endif()
