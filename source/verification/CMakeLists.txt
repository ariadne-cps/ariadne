set(LIBRARY_NAME ariadne-verification)

add_library(${LIBRARY_NAME} OBJECT
    reach_avoid.cpp
    reachability_graph.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    include(CodeCoverage)
    append_coverage_compiler_flags()
endif()
