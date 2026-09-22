set(Python3_FIND_VIRTUALENV FIRST)
set(Python3_FIND_STRATEGY LOCATION)
set(PYBIND11_FINDPYTHON ON)

find_package(Python3 3.11 COMPONENTS Interpreter Development.Module)

option(PYTHON_ENFORCE_TESTS "Enforce Python tests" OFF)

if(Python3_FOUND)

    message(STATUS "Python executable: ${Python3_EXECUTABLE}")
    message(STATUS "Python version:    ${Python3_VERSION}")
    message(STATUS "Python includes:   ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python SOABI:      ${Python3_SOABI}")

    add_subdirectory(pybind11)
    add_subdirectory(source)

    pybind11_add_module(pyariadne
            $<TARGET_OBJECTS:pyariadne-module-obj>
            $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
            $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
            $<TARGET_OBJECTS:pyariadne-extra-submodule-obj>
            )
    target_link_libraries(pyariadne PUBLIC ariadne)

    pybind11_add_module(pyariadne-numeric
            $<TARGET_OBJECTS:pyariadne-numeric-module-obj>
            $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
            $<TARGET_OBJECTS:ariadne-utility>
            $<TARGET_OBJECTS:ariadne-foundations>
            $<TARGET_OBJECTS:ariadne-numeric>
            )
    set_target_properties(pyariadne-numeric PROPERTIES OUTPUT_NAME "pyariadne_numeric")
    target_link_libraries(pyariadne-numeric PUBLIC ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

    pybind11_add_module(pyariadne-algebra
            $<TARGET_OBJECTS:pyariadne-algebra-module-obj>
            $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
            $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
            )
    set_target_properties(pyariadne-algebra PROPERTIES OUTPUT_NAME "pyariadne_algebra")
    target_link_libraries(pyariadne-algebra PUBLIC ariadne)

    set(PYTHON_BINDINGS_INSTALL_DIR "${Python3_SITEARCH}"
            CACHE STRING "Install destination for Python targets")
    install(TARGETS pyariadne DESTINATION "${PYTHON_BINDINGS_INSTALL_DIR}")

    add_subdirectory(tests)

else()
    if(PYTHON_ENFORCE_TESTS)
        message(FATAL_ERROR "Python >= 3.11 with Interpreter and Development.Module components is required when PYTHON_ENFORCE_TESTS=ON.")
    else()
        message(STATUS "Compatible Python development environment not found: Python bindings will not be built.")
    endif()
endif()
