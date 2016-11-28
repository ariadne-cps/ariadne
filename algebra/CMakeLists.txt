add_library(ariadne-algebra OBJECT
    vector.cc
    matrix.cc
    series.cc
    multi_index.cc
    expansion.cc
    differential.cc
    algebra.cc
)

cotire(ariadne-algebra)
