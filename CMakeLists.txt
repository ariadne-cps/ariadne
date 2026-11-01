cmake_minimum_required(VERSION 3.15)
project(Ariadne VERSION 2.5.3)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings when Python is available" ON)

# Configuration is deliberately direct: it bootstraps CMake before Reasoning.
find_package(Git REQUIRED)
if(NOT TARGET ariadne-configuration)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}"
                submodule update --init -- ss/cfg
        RESULT_VARIABLE _configuration_result)
    if(_configuration_result)
        message(FATAL_ERROR "Cannot initialize Configuration submodule")
    endif()
    include("${PROJECT_SOURCE_DIR}/ss/cfg/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration CHAINS
        "configuration" "reasoning|configuration")
endif()

setup_project()
add_subdirectory(ss)

# When included from a parent project, Reasoning exports the 'reasoning'
# interface and its full transitive object list, but not its standalone
# ariadne-reasoning binary. Reuse those same objects for the Ariadne binary.
get_target_property(_ariadne_objects reasoning ARIADNE_PROJECT_OBJECT_TARGETS)
if(NOT _ariadne_objects OR _ariadne_objects STREQUAL "_ariadne_objects-NOTFOUND")
    message(FATAL_ERROR "Reasoning did not export any library object targets")
endif()

add_library(ariadne ${LIBRARY_KIND})
foreach(_object IN LISTS _ariadne_objects)
    target_sources(ariadne PRIVATE $<TARGET_OBJECTS:${_object}>)
endforeach()
target_link_libraries(ariadne PUBLIC reasoning)
if(UNIX AND NOT APPLE)
    target_link_libraries(ariadne PRIVATE dl)
endif()
find_package(Threads REQUIRED)
target_link_libraries(ariadne PRIVATE Threads::Threads)

install(TARGETS ariadne
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)
ariadne_install_dependency_bundle(TARGET reasoning DESTINATION include/ariadne)

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
