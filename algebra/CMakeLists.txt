message(${PROJECT_SOURCE_DIR})

add_library(ariadne-algebra OBJECT
    vector.cc
    matrix.cc
    differential.cc
    algebra.cc
)
