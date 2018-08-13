add_library(ariadne-symbolic OBJECT
    expression.cpp
    expression_set.cpp
)

cotire(ariadne-symbolic)
