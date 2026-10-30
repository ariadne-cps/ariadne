# System definitions, independent of evolution and reachability analysis.
# Build as objects so the staged aggregate 'ariadne' retains its existing
# linking semantics, including on Windows where the aggregate is static.
add_library(ariadne-system OBJECT
    system_function.cpp
    hybrid_constraint_set.cpp
    hybrid_bounded_constraint_set.cpp
    vector_field.cpp
    iterated_map.cpp
    differential_inclusion.cpp
    discrete_location.cpp
    hybrid_automaton.cpp
    hybrid_automaton-composite.cpp
)
set_target_properties(ariadne-system PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-system PUBLIC symbolic expression_set)
