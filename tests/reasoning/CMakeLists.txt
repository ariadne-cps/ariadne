# reasoning layer behavioral tests (moved without changing test logic).
set(UNIT_TESTS
    test_hybrid_reachability
    test_storage_orbit_header
    test_storage_header
    test_enclosure_storage_header
)

# Former dynamics tests for reachability analysis and safety.
set(MIGRATED_REASONING_TESTS
    test_finite_time_reachability
    test_infinite_time_reachability
    test_verify_safety
)
list(APPEND UNIT_TESTS ${MIGRATED_REASONING_TESTS})

foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})
    if("${TEST}" IN_LIST MIGRATED_REASONING_TESTS)
        set_property(TEST ${TEST} PROPERTY LABELS "reasoning")
    else()
        set_property(TEST ${TEST} PROPERTY LABELS "hybrid")
    endif()
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})

if(WIN32)
    set_tests_properties(test_hybrid_reachability PROPERTIES TIMEOUT 240)
    set_tests_properties(test_infinite_time_reachability PROPERTIES TIMEOUT 240)
endif()
