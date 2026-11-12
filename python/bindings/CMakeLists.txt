add_library(pyariadne-ariadne INTERFACE)
target_link_libraries(pyariadne-ariadne INTERFACE pyariadne-reasoning)
ariadne_add_python_module(NAME ariadne MODULE_SOURCE ariadne_module.cpp)
