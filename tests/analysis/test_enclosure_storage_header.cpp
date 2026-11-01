// Verify the analysis-only enclosure-to-storage API compiles as a first include.
#include "analysis/enclosure_storage.hpp"

using namespace Ariadne;

using OuterApproximationFn = Storage (*)(Enclosure const&, Grid const&, Nat);
using AdjoinApproximationFn = Void (*)(Enclosure const&, Storage&, Nat);
using HybridAdjoinApproximationFn = Void (*)(HybridEnclosure const&, HybridStorage&, Nat);

int main() {
    OuterApproximationFn outer = &outer_approximation;
    AdjoinApproximationFn adjoin = &adjoin_outer_approximation_to;
    HybridAdjoinApproximationFn hybrid_adjoin = &adjoin_outer_approximation_to;
    return (outer != nullptr && adjoin != nullptr && hybrid_adjoin != nullptr) ? 0 : 1;
}
