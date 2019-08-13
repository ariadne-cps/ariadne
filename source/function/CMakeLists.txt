add_library(ariadne-function OBJECT
    affine.cpp
    polynomial.cpp
    chebyshev_polynomial.cpp
    formula.cpp
    procedure.cpp
    function.cpp
    function_model.cpp
    affine_model.cpp
    taylor_model.cpp
    taylor_function.cpp
)

cotire(ariadne-function CONFIGURATIONS None)

if(COVERAGE)
    target_link_libraries(ariadne-function PUBLIC coverage_config)
endif()