# analysis layer behavioral tests (moved without changing test logic).
set(UNIT_TESTS
    test_hybrid_reachability
)

foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})
    set_property(TEST ${TEST} PROPERTY LABELS "hybrid")
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})

if(WIN32)
    set_tests_properties(test_hybrid_reachability PROPERTIES TIMEOUT 240)
endif()
