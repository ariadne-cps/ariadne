add_library(ariadne-geometry OBJECT
    point.cpp
    interval.cpp
    box.cpp
    curve.cpp
    function_set.cpp
    binary_tree.cpp
    grid.cpp
    grid_set.cpp
    affine_set.cpp
    enclosure.cpp
    paver.cpp
)

cotire(ariadne-geometry)
