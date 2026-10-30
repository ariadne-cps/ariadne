# Hybrid sets and spaces are an independent representation layer above labelling.
add_library(ariadne-hybrid OBJECT
    discrete_location.cpp
    hybrid_constraint_set.cpp
    hybrid_bounded_constraint_set.cpp
    hybrid_set.cpp
    hybrid_set_adapters.cpp
    hybrid_paving.cpp
    hybrid_graphics.cpp
)
set_target_properties(ariadne-hybrid PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-hybrid PUBLIC labelling)
