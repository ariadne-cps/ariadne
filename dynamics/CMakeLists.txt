message(${PROJECT_SOURCE_DIR})

add_library(ariadne-dynamics OBJECT
    orbit.cc
    map_evolver.cc
    vector_field_evolver.cc
)

cotire(ariadne-dynamics)
