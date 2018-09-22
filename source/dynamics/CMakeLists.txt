add_library(ariadne-dynamics OBJECT
    orbit.cpp
    differential_inclusion.cpp
    map_evolver.cpp
    vector_field_evolver.cpp
    reachability_analyser.cpp
)

cotire(ariadne-dynamics CONFIGURATIONS None)
