set(THREEWATERTANKS_EXAMPLES
    threewatertanks
)

foreach(example ${THREEWATERTANKS_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_dependencies(examples ${THREEWATERTANKS_EXAMPLES})
