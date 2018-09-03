set(UNIT_TESTS
    test_hybrid_automaton
    test_hybrid_evolution    
    test_hybrid_evolver
    test_hybrid_graphics
    test_hybrid_reachability_analyser
    test_hybrid_set
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    add_test(${TEST} ${TEST})
    target_link_libraries(${TEST} ariadne)    
endforeach()

add_dependencies(tests ${UNIT_TESTS})