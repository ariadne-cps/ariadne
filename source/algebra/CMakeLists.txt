add_library(ariadne-algebra OBJECT
    vector.cpp
    matrix.cpp
    series.cpp
    multi_index.cpp
    expansion.cpp
    differential.cpp
    algebra.cpp
)

cotire(ariadne-algebra CONFIGURATIONS None)

if(COVERAGE)
    target_compile_options(ariadne-algebra PUBLIC ${COVERAGE_COMPILER_FLAGS})
endif()