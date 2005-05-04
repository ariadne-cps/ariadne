#include <iostream>

extern const char* filename;

void
test_assert(bool condition, std::string error) {
    if(!condition) {
	std::cout << filename << ": FAILED " << error << '\n';
	exit(1);
    }
}

