/***************************************************************************
 *            solvers/solver_attributes.hpp
 *
 *  Copyright  2011-26  Pieter Collins
 *
 ****************************************************************************/

#ifndef ARIADNE_SOLVER_ATTRIBUTES_HPP
#define ARIADNE_SOLVER_ATTRIBUTES_HPP

#include "utility/attribute.hpp"
#include "numeric/builtin.hpp"

namespace Ariadne {


struct MaximumError : Attribute<ApproximateDouble> {
    MaximumError(ApproximateDouble x) : Attribute<ApproximateDouble>(x) { }
    MaximumError(double x) : Attribute<ApproximateDouble>(x) { }
};
struct SweepThreshold : Attribute<ApproximateDouble> { using Attribute<ApproximateDouble>::Attribute; };
struct MaximumNumericTypeOfSteps : Attribute<Nat> { using Attribute<Nat>::Attribute; };

static const Generator<MaximumError> maximum_error = Generator<MaximumError>();
static const Generator<SweepThreshold> sweep_threshold = Generator<SweepThreshold>();
static const Generator<MaximumNumericTypeOfSteps> maximum_number_of_steps = Generator<MaximumNumericTypeOfSteps>();

} // namespace Ariadne

#endif /* ARIADNE_SOLVER_ATTRIBUTES_HPP */
