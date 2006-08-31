#include <iostream>

void
test_assert(bool condition, std::string error) {
  if(!condition) {
    std::cout << "FAILED " << error << '\n' << std::flush;
    exit(1);
  }
}

void test_assert(bool condition) {
  if(!condition) {
    std::cout << "FAILED" << '\n' << std::flush;
    exit(1);
  }
}
