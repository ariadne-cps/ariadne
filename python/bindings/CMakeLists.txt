include_directories(${Python3_INCLUDE_DIRS})
get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-bindings-obj OBJECT
    reasoning_submodule.cpp
)
add_library(pyariadne-module-obj OBJECT ariadne_module.cpp)
foreach(target pyariadne-bindings-obj pyariadne-module-obj)
    target_include_directories(${target} PRIVATE
        "${PYTHON_COMMON_SOURCE_DIR}/include"
        $<TARGET_PROPERTY:pyariadne-algebra,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(${target} PRIVATE ariadne pyariadne-evolution pyariadne-python-common pybind11::module)
    set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()
