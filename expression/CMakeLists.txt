add_library(ariadne-expression OBJECT
    formula.cc
    expression.cc
    expression_set.cc
)

cotire(ariadne-expression)
