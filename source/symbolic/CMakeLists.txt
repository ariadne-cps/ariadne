add_library(ariadne-symbolic OBJECT
    space.cpp
    expression.cpp
    expression_set.cpp
)

cotire(ariadne-symbolic CONFIGURATIONS None)

if(COVERAGE)
    target_link_libraries(ariadne-symbolic PUBLIC coverage_config)
endif()
