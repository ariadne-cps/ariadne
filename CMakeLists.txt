cmake_minimum_required(VERSION 3.15)
project(Ariadne VERSION 3.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings when Python is available" ON)

# Build must be bootstrapped before Reasoning so the CMake helpers are available.
find_package(Git REQUIRED)
if(NOT TARGET ariadne-build)
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env "LC_ALL=C"
                "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}"
                submodule update --init -- ss/bui
        RESULT_VARIABLE _build_result)
    if(_build_result)
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

# Use Build's standard library assembly and installation logic.
get_target_property(_ariadne_objects reasoning ARIADNE_PROJECT_OBJECT_TARGETS)
if(NOT _ariadne_objects)
    message(FATAL_ERROR "Reasoning did not export any library object targets")
endif()

setup_project_library(
    TARGET ariadne
    OBJECTS ${_ariadne_objects}
    LINK_LIBRARIES reasoning)
# setup_project_library creates the archive only for a root project.
# Keep the previous behavior when Ariadne is included by another project.
if(NOT TARGET ariadne)
    add_library(ariadne ${LIBRARY_KIND})
    foreach(_object IN LISTS _ariadne_objects)
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
ariadne_register_public_headers(
    TARGET ariadne
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne.hpp")

if(UNIX AND NOT APPLE)
    target_link_libraries(ariadne PRIVATE dl)
endif()
find_package(Threads REQUIRED)
target_link_libraries(ariadne PRIVATE Threads::Threads)

ariadne_install_dependency_bundle(TARGET ariadne DESTINATION include/ariadne)

# Reasoning registers every upstream Python module, including SAT and SMT.
# Compose one extension from the registered entry points and the objects
# exported by pyariadne-reasoning. No local binding sources are needed.
if(BUILD_PYTHON_BINDINGS AND TARGET pyariadne-reasoning)
    find_package(Python3 3.11 COMPONENTS Interpreter Development.Module QUIET)
    if(Python3_FOUND)
        get_property(_modules GLOBAL PROPERTY ARIADNE_PYTHON_MODULES)
        if(NOT "reasoning" IN_LIST _modules)
            message(FATAL_ERROR "Reasoning Python module was not registered")
        endif()
        get_target_property(_python_common_source pyariadne-python-common SOURCE_DIR)
        set(_entry "#include \"pybind11.hpp\"\n\n")
        foreach(_module IN LISTS _modules)
            string(APPEND _entry "void ${_module}_module(pybind11::module&);\n")
        endforeach()
        string(APPEND _entry "\nPYBIND11_MODULE(pyariadne, module) {\n")
        foreach(_module IN LISTS _modules)
            string(APPEND _entry "    ${_module}_module(module);\n")
        endforeach()
        string(APPEND _entry "}\n")
        file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pyariadne_entry.cpp"
            CONTENT "${_entry}")
        pybind11_add_module(pyariadne
            "${CMAKE_CURRENT_BINARY_DIR}/pyariadne_entry.cpp")
        target_include_directories(pyariadne PRIVATE
            "${_python_common_source}/include")
        foreach(_object IN LISTS _ariadne_objects)
            target_sources(pyariadne PRIVATE $<TARGET_OBJECTS:${_object}>)
        endforeach()
        target_link_libraries(pyariadne PRIVATE ariadne pyariadne-reasoning)
        if(WIN32)
            set(_python_install_dir "Lib/site-packages")
        else()
            set(_python_install_dir
                "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
        endif()
        install(TARGETS pyariadne DESTINATION "${_python_install_dir}")
        add_dependencies(bindings pyariadne)
    endif()
endif()
