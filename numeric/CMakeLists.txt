message(${PROJECT_SOURCE_DIR})

add_library(ariadne-numeric OBJECT
    logical.cc
    integer.cc
    rational.cc
    decimal.cc
    real.cc
    float.cc
    float-approximate.cc
    float-validated.cc
    floatmp.cc
)
