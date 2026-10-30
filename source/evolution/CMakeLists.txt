# Evolution of continuous and hybrid systems, including orbit storage.
add_library(ariadne-evolution OBJECT
    ../dynamics/enclosure.cpp
    ../dynamics/orbit.cpp
    ../dynamics/differential_inclusion_evolver.cpp
    ../dynamics/differential_inclusion_grid.cpp
    ../dynamics/iterated_map_evolver.cpp
    ../dynamics/vector_field_evolver.cpp
    ../dynamics/vector_field_simulator.cpp
    ../hybrid/hybrid_set_orbit.cpp
    ../hybrid/hybrid_enclosure.cpp
    ../hybrid/hybrid_simulator.cpp
    ../hybrid/hybrid_evolver.cpp
    hybrid_paving_evolution.cpp
)
set_target_properties(ariadne-evolution PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-evolution PUBLIC ariadne-system integration equation_solving)
