set(WATERTANK_EXAMPLES 
    watertank
    watertank-monolithic
)

foreach(example ${WATERTANK_EXAMPLES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()