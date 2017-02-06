add_library(ariadne-numeric OBJECT
    dyadic.cpp
    logical.cpp
    operators.cpp
    integer.cpp
    rational.cpp
    decimal.cpp
    number.cpp
    real.cpp
    float64.cpp
    floatmp.cpp
    float-user.cpp
)

cotire(ariadne-numeric)
