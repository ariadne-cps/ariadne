add_library(ariadne-hybrid OBJECT
    discrete_location.cc
    hybrid_set.cc
    hybrid_enclosure.cc
    hybrid_automaton.cc
    hybrid_automaton-composite.cc
    hybrid_automaton-restrictive.cc
    hybrid_simulator.cc
    hybrid_evolver.cc
    hybrid_reachability_analyser.cc
    hybrid_graphics.cc
)

cotire(ariadne-hybrid)
