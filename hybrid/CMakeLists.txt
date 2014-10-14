message(${PROJECT_SOURCE_DIR})

add_library(ariadne-hybrid OBJECT
    hybrid_set.cc
    hybrid_enclosure.cc
    hybrid_automaton-monolithic.cc
    hybrid_automaton-composite.cc
    hybrid_automaton-restrictive.cc
    hybrid_simulator.cc
    hybrid_evolver.cc
    hybrid_reachability_analyser.cc
    hybrid_graphics.cc
)
