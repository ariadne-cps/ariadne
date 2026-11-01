// Regression: analysis storage must be self-contained when first included.
#include "reasoning/storage.hpp"
#include <type_traits>

static_assert(std::is_same<Ariadne::RealSpace, Ariadne::Space<Ariadne::Real>>::value,
              "Analysis storage must see the complete symbolic state space");

int main() { return 0; }
