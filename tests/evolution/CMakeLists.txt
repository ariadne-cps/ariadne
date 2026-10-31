# evolution layer behavioral tests (moved without changing test logic).
set(UNIT_TESTS
    test_hybrid_enclosure
    test_hybrid_evolution
    test_hybrid_evolver
    test_hybrid_simulator
    test_first_order_pde
    test_first_order_pde_evolution_header
    test_orbit_evolution_header
)

foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})
    set_property(TEST ${TEST} PROPERTY LABELS "hybrid")
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})

# A source-level guard, independent of the compiled numerical tests.
add_test(NAME test_evolution_dependency_boundary
    COMMAND ${CMAKE_COMMAND}
        "-DEVOLUTION_SOURCE_DIR=${PROJECT_SOURCE_DIR}/source/evolution"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_dependency_boundary.cmake")
