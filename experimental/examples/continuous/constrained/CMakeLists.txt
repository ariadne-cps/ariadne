set(EXAMPLES_CONSTRAINED
    constrained_vanderpol
)

foreach(example ${EXAMPLES_CONSTRAINED})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()
