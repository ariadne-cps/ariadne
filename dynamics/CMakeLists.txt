add_library(ariadne-dynamics OBJECT
    orbit.cpp
    map_evolver.cpp
    vector_field_evolver.cpp
)

cotire(ariadne-dynamics)
