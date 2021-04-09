set(LIBRARY_NAME ariadne-dynamics)

add_library(${LIBRARY_NAME} OBJECT
    enclosure.cpp
    orbit.cpp
    differential_inclusion.cpp
    differential_inclusion_evolver.cpp
    first_order_pde.cpp
    iterated_map.cpp
    iterated_map_evolver.cpp
    vector_field.cpp
    vector_field_evolver.cpp
    vector_field_simulator.cpp
    reachability_analyser.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
