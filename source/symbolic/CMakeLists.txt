add_library(ariadne-symbolic OBJECT
    space.cpp
    expression.cpp
    expression_set.cpp
)

cotire(ariadne-symbolic CONFIGURATIONS None)

if(COVERAGE)
    target_compile_options(ariadne-symbolic PUBLIC ${COVERAGE_COMPILER_FLAGS})
endif()
