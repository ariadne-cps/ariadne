add_library(ariadne-dynamics OBJECT
    enclosure.cpp
    orbit.cpp
    inclusion_vector_field.cpp
    inclusion_evolver.cpp
    first_order_pde.cpp
    map_evolver.cpp
    vector_field_evolver.cpp
    reachability_analyser.cpp
)

cotire(ariadne-dynamics CONFIGURATIONS None)

if(COVERAGE)
    target_link_libraries(ariadne-dynamics PUBLIC coverage_config)
endif()
