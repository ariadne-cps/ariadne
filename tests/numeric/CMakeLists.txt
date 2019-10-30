set(UNIT_TESTS
    test_logical
    test_integer
    test_dyadic
    test_rational
    test_real
    test_number
    test_rounding
    test_float
    test_float_approximation
    test_float_lower_upper_bound
    test_float_bounds
    test_float_ball
    test_float_value
    test_complex
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    if(COVERAGE)
        target_compile_options(${TEST} PUBLIC ${COVERAGE_COMPILER_FLAGS})
    endif()
    target_link_libraries(${TEST} ariadne-core)
    add_test(${TEST} ${TEST})
    set_property(TEST ${TEST} PROPERTY LABELS "numeric;primary;core;kernel")
endforeach()

add_dependencies(tests-core ${UNIT_TESTS})
add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})

#add_library(check_numeric OBJECT check_numeric.cpp)
#add_dependencies(tests-core check_numeric)
#add_dependencies(tests-kernel check_numeric)
#add_dependencies(tests check_numeric)
