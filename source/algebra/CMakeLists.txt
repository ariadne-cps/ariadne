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
    target_link_libraries(ariadne-algebra PUBLIC coverage_config)
endif()