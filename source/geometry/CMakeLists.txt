add_library(ariadne-geometry OBJECT
    point.cpp
    interval.cpp
    box.cpp
    curve.cpp
    function_set.cpp
    binary_tree.cpp
    grid.cpp
    grid_paving.cpp
    paver.cpp
    affine_set.cpp
)

cotire(ariadne-geometry CONFIGURATIONS None)

if(COVERAGE)
    target_link_libraries(ariadne-geometry PUBLIC coverage_config)
endif()
