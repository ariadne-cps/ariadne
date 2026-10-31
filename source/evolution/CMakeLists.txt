# Evolution of continuous and hybrid systems, including orbit storage.
add_library(ariadne-evolution OBJECT
    enclosure.cpp
    orbit.cpp
    differential_inclusion_evolver.cpp
    iterated_map_evolver.cpp
    vector_field_evolver.cpp
    vector_field_simulator.cpp
    first_order_pde.cpp
    hybrid_orbit.cpp
    hybrid_enclosure.cpp
    hybrid_simulator.cpp
    hybrid_evolver.cpp
)
set_target_properties(ariadne-evolution PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-evolution PUBLIC system integration equation_solving)
