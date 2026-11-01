# Reachability and higher-level reasoning, never referenced by lower layers.
add_library(ariadne-reasoning OBJECT
    reachability_analyser.cpp
    hybrid_reachability_analyser.cpp
    hybrid_storage_orbit.cpp
    storage_orbit.cpp
    enclosure_storage.cpp
)
set_target_properties(ariadne-reasoning PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(ariadne-reasoning PUBLIC evolution smt sat)
