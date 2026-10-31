# Reachability and higher-level analysis, never referenced by lower layers.
add_library(ariadne-analysis OBJECT
    ../dynamics/reachability_analyser.cpp
    ../dynamics/first_order_pde.cpp
    hybrid_reachability_analyser.cpp
    hybrid_storage_orbit.cpp
)
set_target_properties(ariadne-analysis PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-analysis PUBLIC ariadne-evolution smt sat)
