# Reasoning has already registered all transitive Python modules.
set(Python3_FIND_VIRTUALENV FIRST)
set(Python3_FIND_STRATEGY LOCATION)
set(PYBIND11_FINDPYTHON ON)
find_package(Python3 3.11 COMPONENTS Interpreter Development.Module REQUIRED)

add_subdirectory(bindings)

# Build's finalizer gathers the complete C++ object set from this target.
get_target_property(_objects reasoning ARIADNE_PROJECT_OBJECT_TARGETS)
set_property(TARGET ariadne PROPERTY ARIADNE_PROJECT_OBJECT_TARGETS "${_objects}")
ariadne_finalize_python_module(NAME ariadne)

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
