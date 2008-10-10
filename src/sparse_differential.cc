#include "sparse_differential.h"


#include "differential_vector.h"

namespace Ariadne {
template class SparseDifferential<Float>;
template class DifferentialVector< SparseDifferential<Float> >;
}
