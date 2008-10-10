#include "dense_differential.h"

#include "differential_vector.h"

namespace Ariadne {
template class DenseDifferential<Float>;
template class DifferentialVector< DenseDifferential<Float> >;
}
