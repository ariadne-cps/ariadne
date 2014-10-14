message(${PROJECT_SOURCE_DIR})

add_library(ariadne-expression OBJECT
    operators.cc
    expression.cc
    expression_set.cc
)
