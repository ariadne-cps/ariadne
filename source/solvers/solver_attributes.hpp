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

using Utility::Attribute;
using Utility::Generator;

struct MaximumError : Attribute<ApproximateDouble> {
    MaximumError(ApproximateDouble x) : Attribute<ApproximateDouble>(x) { }
    MaximumError(double x) : Attribute<ApproximateDouble>(x) { }
};
struct SweepThreshold : Attribute<ApproximateDouble> { using Attribute<ApproximateDouble>::Attribute; };
struct MaximumNumericTypeOfSteps : Attribute<Nat> { using Attribute<Nat>::Attribute; };

static const Generator<MaximumError> maximum_error = Generator<MaximumError>();
static const Generator<SweepThreshold> sweep_threshold = Generator<SweepThreshold>();
static const Generator<MaximumNumericTypeOfSteps> maximum_number_of_steps = Generator<MaximumNumericTypeOfSteps>();

struct Capacity : Attribute<SizeType> { };
static const Generator<Capacity> capacity = Generator<Capacity>();
struct Size : Attribute<SizeType> { };
static const Generator<Size> size = Generator<Size>();
struct ResultSize : Attribute<SizeType> { };
static const Generator<ResultSize> result_size = Generator<ResultSize>();
struct ArgumentSize : Attribute<SizeType> { };
static const Generator<ArgumentSize> argument_size = Generator<ArgumentSize>();
struct Degree : Attribute<DegreeType> { };
static const Generator<Degree> degree = Generator<Degree>();

} // namespace Ariadne

#endif /* ARIADNE_SOLVER_ATTRIBUTES_HPP */
