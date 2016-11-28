add_library(ariadne-numeric OBJECT
    dyadic.cc
    logical.cc
    operators.cc
    integer.cc
    rational.cc
    decimal.cc
    number.cc
    real.cc
    float64.cc
    floatmp.cc
    float-user.cc
)

cotire(ariadne-numeric)
