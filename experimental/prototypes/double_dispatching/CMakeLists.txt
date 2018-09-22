cmake_minimum_required (VERSION 3.9)

set(CMAKE_CXX_STANDARD 17)

add_executable(double_dispatching double_dispatching.cpp)
target_link_libraries(double_dispatching ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})