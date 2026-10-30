# System specifications and dynamics, layered on hybrid set representations.
add_library(ariadne-system OBJECT
    system_function.cpp
    vector_field.cpp
    iterated_map.cpp
    differential_inclusion.cpp
    hybrid_automaton.cpp
    hybrid_automaton-composite.cpp
)
set_target_properties(ariadne-system PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-system PUBLIC ariadne-hybrid)
