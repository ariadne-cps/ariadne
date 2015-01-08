message(${PROJECT_SOURCE_DIR})

add_library(ariadne-algebra OBJECT
    vector.cc
    matrix.cc
    multi_index.cc
    expansion.cc
    differential.cc
    algebra.cc
)
