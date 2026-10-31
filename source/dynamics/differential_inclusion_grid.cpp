// The default spatial grid belongs to paving, not to system definitions.
// Keep the legacy DifferentialInclusion::grid() API in the dynamics adapter
// without making ariadne-system depend on the paving repository.
#include "dynamics/differential_inclusion.hpp"

namespace Ariadne {

Grid DifferentialInclusion::grid() const {
    return Grid(this->dimension());
}

} // namespace Ariadne
