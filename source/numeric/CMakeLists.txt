set(LIBRARY_NAME ariadne-numeric)

add_library(${LIBRARY_NAME} OBJECT
    operators.cpp
    int.cpp
    double.cpp
    integer.cpp
    dyadic.cpp
    decimal.cpp
    rational.cpp
    real.cpp
    complex.cpp
    number.cpp
    floatdp.cpp
    floatmp.cpp
    float_approximation.cpp
    float_lower_bound.cpp
    float_upper_bound.cpp
    float_bounds.cpp
    float_ball.cpp
    float_error.cpp
    float_literals.cpp

    concepts.cpp
    logical_sequence.cpp
)

target_link_libraries(${LIBRARY_NAME}
    PUBLIC
        foundation
        utility
        GMP::GMP
        MPFR::MPFR
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR NOT APPLE)
    target_compile_options(${LIBRARY_NAME} PRIVATE -frounding-math)
endif()
