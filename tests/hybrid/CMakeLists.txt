set(UNIT_TESTS
    test_hybrid_automaton
    test_hybrid_evolution    
    test_hybrid_evolver
    test_hybrid_graphics
    test_hybrid_set
    test_hybrid_simulator
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    if(COVERAGE)
        target_compile_options(${TEST} PUBLIC ${COVERAGE_COMPILER_FLAGS})
    endif()    
    target_link_libraries(${TEST} ariadne)  
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "hybrid;primary")           
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})