set(UNIT_TESTS
    test_crossings_issue
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne)       
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "hybrid;secondary") 
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})