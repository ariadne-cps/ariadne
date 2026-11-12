# Reasoning already registers all transitively available Python modules.
find_package(Python3 3.11 COMPONENTS Interpreter Development.Module REQUIRED)
add_subdirectory(bindings)

if(WIN32)
    set(_python_install_dir "Lib/site-packages")
else()
    set(_python_install_dir
        "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
endif()
set(PYTHON_BINDINGS_INSTALL_DIR "${_python_install_dir}"
    CACHE STRING "Python extension install directory relative to prefix")
install(TARGETS pyariadne DESTINATION "${PYTHON_BINDINGS_INSTALL_DIR}")
add_dependencies(bindings pyariadne)
