// Regression: the system-level PDE model header is independently includable.
// In particular, algebra/matrix.hpp needs a complete Rational definition.
#include "system/first_order_pde.hpp"

int main() {
    return 0;
}
