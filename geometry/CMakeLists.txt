message(${PROJECT_SOURCE_DIR})

add_library(ariadne-geometry OBJECT
    point.cc
    interval.cc
    box.cc
    curve.cc
    function_set.cc
    grid_set.cc
    affine_set.cc
    enclosure.cc
    paver.cc
)
