set(UNIT_TESTS
    test_dyadic
    test_float
    test_integer
    test_logical
    test_number
    test_rational
    test_real    
    test_rounding
    test_validated_float
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-core)    
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "numeric;primary;core;kernel")      
endforeach()

add_dependencies(tests-core ${UNIT_TESTS})
add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})