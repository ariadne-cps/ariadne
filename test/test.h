#include <iostream>

void
test_assert(bool condition, std::string error) {
  if(!condition) {
    std::cerr << "FAILED " << error << '\n' << std::flush;
    exit(1);
  }
}

void test_assert(bool condition) {
  if(!condition) {
    std::cerr << "FAILED" << '\n' << std::flush;
    exit(1);
  }
}
