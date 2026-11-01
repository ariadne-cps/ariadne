// The grid-storage Orbit specialization is owned by analysis.
#include "reasoning/storage_orbit.hpp"
#include <type_traits>

static_assert(std::is_same<Ariadne::Orbit<Ariadne::Storage>::EnclosureListType,
                           Ariadne::Storage>::value,
              "Orbit<Storage> must be complete in the analysis header");

int main() {
    return 0;
}
