set(UNIT_TESTS
    test_continuous_evolution
    test_inclusion_vector_field
    test_inclusion_evolver
    test_discrete_evolution
    test_reachability_analyser
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    if(COVERAGE)
        target_compile_options(${TEST} PUBLIC ${COVERAGE_COMPILER_FLAGS})
    endif()    
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "dynamics;primary")
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
