get_property(_modules GLOBAL PROPERTY ARIADNE_PYTHON_MODULES)
if(NOT "reasoning" IN_LIST _modules)
    message(FATAL_ERROR "Reasoning Python module was not registered")
endif()

foreach(_module IN LISTS _modules)
    string(APPEND ARIADNE_PY_DECLARATIONS
        "void ${_module}_module(pybind11::module&);\n")
    string(APPEND ARIADNE_PY_INITIALIZERS
        "    ${_module}_module(module);\n")
endforeach()
configure_file(pyariadne_entry.cpp.in pyariadne_entry.cpp @ONLY)

pybind11_add_module(pyariadne "${CMAKE_CURRENT_BINARY_DIR}/pyariadne_entry.cpp")
target_link_libraries(pyariadne PRIVATE
    ariadne pyariadne-reasoning pyariadne-python-common)
get_target_property(_objects reasoning ARIADNE_PROJECT_OBJECT_TARGETS)
foreach(_object IN LISTS _objects)
    target_sources(pyariadne PRIVATE $<TARGET_OBJECTS:${_object}>)
endforeach()
