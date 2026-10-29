set(LIBRARY_NAME ariadne-hybrid)

# Transitional object library: files currently mix system/evolution/analysis.
add_library(${LIBRARY_NAME} OBJECT
    discrete_location.cpp
    hybrid_set.cpp
    hybrid_paving.cpp
    hybrid_enclosure.cpp
    hybrid_automaton.cpp
    hybrid_automaton-composite.cpp
    hybrid_simulator.cpp
    hybrid_evolver.cpp
    hybrid_reachability_analyser.cpp
    hybrid_graphics.cpp
)
target_link_libraries(${LIBRARY_NAME} PRIVATE ariadne-analysis)
