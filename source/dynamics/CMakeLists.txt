set(LIBRARY_NAME ariadne-dynamics)

# Transitional object library: remaining files are evolution/analysis.
add_library(${LIBRARY_NAME} OBJECT
    enclosure.cpp
    orbit.cpp
    differential_inclusion_evolver.cpp
    differential_inclusion_grid.cpp
    first_order_pde.cpp
    iterated_map_evolver.cpp
    vector_field_evolver.cpp
    vector_field_simulator.cpp
    reachability_analyser.cpp
)
target_link_libraries(${LIBRARY_NAME} PRIVATE ariadne-analysis)
