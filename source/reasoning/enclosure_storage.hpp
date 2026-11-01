// Bridges evolution enclosures to reasoning-owned grid storage.
// This interface must not be declared in evolution/enclosure.hpp.
#ifndef ARIADNE_REASONING_ENCLOSURE_STORAGE_HPP
#define ARIADNE_REASONING_ENCLOSURE_STORAGE_HPP

#include "reasoning/storage.hpp"
#include "evolution/enclosure.hpp"

namespace Ariadne {

// Adjoin the outer approximation of an evolved enclosure to analysis storage.
Void adjoin_outer_approximation_to(Enclosure const& enclosure, Storage& storage, Nat fineness);

// Construct analysis storage containing an outer approximation of the enclosure.
Storage outer_approximation(Enclosure const& enclosure, Grid const& grid, Nat fineness);

class HybridEnclosure;
class HybridStorage;

// Preserve the hybrid paving semantics without exposing reasoning Storage in evolution.
Void adjoin_outer_approximation_to(HybridEnclosure const& enclosure, HybridStorage& storage, Nat fineness);

} // namespace Ariadne

#endif // ARIADNE_REASONING_ENCLOSURE_STORAGE_HPP
