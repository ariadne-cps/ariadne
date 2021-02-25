cmake_minimum_required(VERSION 3.9)

project(PythonBindings VERSION 2.1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(Ariadne REQUIRED)

include_directories(SYSTEM ${GMP_INCLUDES} ${MPFR_INCLUDES} ${ARIADNE_INCLUDE_DIRS})

find_package(PythonLibs REQUIRED)

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/pybind11/CMakeLists.txt)
    message(STATUS "pybind11 dependency not available.")
    find_package(Git)
    if(GIT_FOUND)
        message(STATUS "Downloading pybind11 repository...")
        if (NOT EXISTS ${PROJECT_SOURCE_DIR}/.git) # Manages the case when an archive is used
            execute_process(COMMAND git init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
            execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/pybind11 OUTPUT_QUIET ERROR_QUIET)
            execute_process(COMMAND git submodule add -f https://github.com/pybind/pybind11 pybind11 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
            execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
        else() # When using a clone
            execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
        endif()
        if(NOT EXISTS ${PROJECT_SOURCE_DIR}/pybind11/CMakeLists.txt)
            message(ERROR "pybind11 download failure.")
        else()
            message(STATUS "Downloaded pybind11 dependency successfully.")
        endif()
    else()
        message(ERROR "Git not found, pybind11 dependency could not be downloaded.")
    endif()
endif()

add_subdirectory(pybind11)
add_subdirectory(source)

add_library(pyariadne MODULE
    $<TARGET_OBJECTS:pyariadne-module-obj>
    $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
    $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
    $<TARGET_OBJECTS:pyariadne-extra-submodule-obj>
)
set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne" SUFFIX ".so")
target_link_libraries (pyariadne PUBLIC ${ARIADNE_LIBRARIES} PRIVATE pybind11::module)
