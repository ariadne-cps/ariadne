add_library(ariadne-utility OBJECT
    utility.cpp
)

if(COVERAGE)
    target_compile_options(ariadne-utility PUBLIC ${COVERAGE_COMPILER_FLAGS})
endif()