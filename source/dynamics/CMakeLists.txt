add_library(ariadne-dynamics OBJECT
    orbit.cpp
    inclusion_vector_field.cpp
    differential_inclusion.cpp
    map_evolver.cpp
    vector_field_evolver.cpp
    reachability_analyser.cpp
)

cotire(ariadne-dynamics CONFIGURATIONS None)

if(COVERAGE)
    target_compile_options(ariadne-dynamics PUBLIC ${COVERAGE_COMPILER_FLAGS})
endif()