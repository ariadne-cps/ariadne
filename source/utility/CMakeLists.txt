add_library(ariadne-utility OBJECT
    utility.cpp
)

target_link_libraries(ariadne-utility PUBLIC coverage_config)