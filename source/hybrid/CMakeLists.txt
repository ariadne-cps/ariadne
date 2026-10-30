set(LIBRARY_NAME ariadne-hybrid)

# Transitional object library: hybrid automata and related types migrate in stages.
add_library(${LIBRARY_NAME} OBJECT
    hybrid_set.cpp
    hybrid_paving.cpp
    hybrid_enclosure.cpp
    hybrid_simulator.cpp
    hybrid_evolver.cpp
    hybrid_reachability_analyser.cpp
    hybrid_graphics.cpp
)
target_link_libraries(${LIBRARY_NAME} PRIVATE ariadne-analysis)
