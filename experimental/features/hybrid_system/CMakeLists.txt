add_library(hybrid_system SHARED hybrid_system.cpp)
target_include_directories(hybrid_system PUBLIC ../../../source/hybrid)
target_link_libraries(hybrid_system ariadne)

add_executable(test_hybrid_system test_hybrid_system.cpp)
target_include_directories(test_hybrid_system PUBLIC ../../../tests/hybrid ../../../source/hybrid)
target_link_libraries(test_hybrid_system hybrid_system)