add_library(ariadne-numeric OBJECT
    dyadic.cpp
    logical.cpp
    operators.cpp
    integer.cpp
    rational.cpp
    decimal.cpp
    number.cpp
    real.cpp
    complex.cpp
    floatdp.cpp
    floatmp.cpp
    float_approximation.cpp
    float_lower_bound.cpp
    float_upper_bound.cpp
    float_bounds.cpp
    float_ball.cpp
    float_value.cpp
    float_error.cpp
    float_literals.cpp
)

cotire(ariadne-numeric CONFIGURATIONS None)

if(COVERAGE)
     target_link_libraries(ariadne-numeric PUBLIC coverage_config)
endif()
