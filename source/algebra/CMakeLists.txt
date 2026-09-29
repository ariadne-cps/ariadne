set(LIBRARY_NAME ariadne-algebra)

add_library(${LIBRARY_NAME} OBJECT
    vector.cpp
    matrix.cpp
    series.cpp
    multi_index.cpp
    expansion.cpp
    differential.cpp
    algebra.cpp
)


