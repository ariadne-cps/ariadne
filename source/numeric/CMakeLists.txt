add_library(ariadne-numeric OBJECT
    dyadic.cpp
    logical.cpp
    operators.cpp
    integer.cpp
    rational.cpp
    decimal.cpp
    number.cpp
    real.cpp
    floatdp.cpp
    floatmp.cpp
    float-user.cpp
)

cotire(ariadne-numeric CONFIGURATIONS None)

if(COVERAGE)
    target_compile_options(ariadne-numeric PUBLIC ${COVERAGE_COMPILER_FLAGS})
endif()
