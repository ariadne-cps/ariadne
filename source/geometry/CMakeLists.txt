set(LIBRARY_NAME ariadne-geometry)

add_library(${LIBRARY_NAME} OBJECT
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
    measurable_set.cpp

    geometry.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
