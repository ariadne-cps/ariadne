add_library(ariadne-utility OBJECT
    utility.cpp
)

if(COVERAGE)
    target_link_libraries(ariadne-utility PUBLIC coverage_config)
endif()