add_library(ariadne-hybrid OBJECT
    discrete_location.cpp
    hybrid_set.cpp
    hybrid_paving.cpp
    hybrid_enclosure.cpp
    hybrid_automaton.cpp
    hybrid_automaton-composite.cpp
    hybrid_automaton-restrictive.cpp
    hybrid_simulator.cpp
    hybrid_evolver.cpp
    hybrid_reachability_analyser.cpp
    hybrid_graphics.cpp
)

cotire(ariadne-hybrid CONFIGURATIONS None)

if(COVERAGE)
    target_link_libraries(ariadne-hybrid PUBLIC coverage_config)
endif()
