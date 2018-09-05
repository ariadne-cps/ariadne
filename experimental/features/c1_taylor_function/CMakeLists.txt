add_library(c1_taylor_function SHARED c1_taylor_function.cpp)
target_include_directories(c1_taylor_function PUBLIC ../../../source/function)
target_link_libraries(c1_taylor_function ariadne)

add_executable(test_c1_taylor_function test_c1_taylor_function.cpp)
target_include_directories(test_c1_taylor_function PUBLIC ../../../tests/function ../../../source/function)
target_link_libraries(test_c1_taylor_function c1_taylor_function)