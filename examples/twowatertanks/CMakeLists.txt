set(TWOWATERTANKS_EXAMPLES
    twowatertanks
)

foreach(example ${TWOWATERTANKS_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples ${TWOWATERTANKS_EXAMPLES})
