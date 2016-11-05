message(${PROJECT_SOURCE_DIR})

add_library(ariadne-function OBJECT
    affine.cc
    polynomial.cc
    function.cc
    function_model.cc
    affine_model.cc
    taylor_model.cc
    taylor_function.cc
)
