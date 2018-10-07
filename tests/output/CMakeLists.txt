set(UNIT_TESTS
    test_graphics
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    if(COVERAGE)
        target_compile_options(${TEST} PUBLIC ${COVERAGE_COMPILER_FLAGS})
    endif()    
    add_test(${TEST} ${TEST})
    target_link_libraries(${TEST} ariadne-kernel) 
    set_property(TEST ${TEST} PROPERTY LABELS "output;primary;kernel")         
endforeach()

add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})