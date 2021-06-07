find_package(Python3 COMPONENTS Interpreter Development)

if(PYTHON3_FOUND)
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/python/pybind11/CMakeLists.txt)
        message(STATUS "pybind11 dependency not available.")
        find_package(Git)
        if(GIT_FOUND)
            message(STATUS "Downloading pybind11 repository...")
            if (NOT EXISTS ${PROJECT_SOURCE_DIR}/.git) # Manages the case when an archive is used
                execute_process(COMMAND git init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
                execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/python/pybind11 OUTPUT_QUIET ERROR_QUIET)
                execute_process(COMMAND git submodule add -f https://github.com/pybind/pybind11 python/pybind11 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
                execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
            else() # When using a clone
                execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
            endif()
            if(NOT EXISTS ${PROJECT_SOURCE_DIR}/python/pybind11/CMakeLists.txt)
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
    set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "pyariadne" SUFFIX ".so")
    target_link_libraries (pyariadne PUBLIC ariadne PRIVATE pybind11::module)

    add_library(pyariadne-numeric MODULE
            $<TARGET_OBJECTS:pyariadne-numeric-module-obj>
            $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
            $<TARGET_OBJECTS:ariadne-utility>
            $<TARGET_OBJECTS:ariadne-numeric>
            )
    set_target_properties(pyariadne-numeric PROPERTIES PREFIX "" OUTPUT_NAME "numeric" SUFFIX ".so")
    target_link_libraries (pyariadne-numeric PUBLIC ${MPFR_LIBRARIES} ${GMP_LIBRARIES} PRIVATE pybind11::module)

    add_library(pyariadne-algebra MODULE
            $<TARGET_OBJECTS:pyariadne-algebra-module-obj>
            $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
            $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
            )
    set_target_properties(pyariadne-algebra PROPERTIES PREFIX "" OUTPUT_NAME "algebra" SUFFIX ".so")
    target_link_libraries (pyariadne-algebra PUBLIC ariadne PRIVATE pybind11::module)

    set(PYTHON_BINDINGS_INSTALL_DIR "${Python3_SITEARCH}"
            CACHE STRING "Install destination for Python targets")
    install(TARGETS pyariadne DESTINATION "${PYTHON_BINDINGS_INSTALL_DIR}")

    add_subdirectory(tests)

else()
    message(STATUS "Python3 not found: Python bindings module 'pyariadne' will not be built.")
endif()
