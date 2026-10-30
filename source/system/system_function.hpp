#ifndef ARIADNE_SYSTEM_SYSTEM_FUNCTION_HPP
#define ARIADNE_SYSTEM_SYSTEM_FUNCTION_HPP

#include "function/function.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"

namespace Ariadne {

//! Order algebraic assignments by dependencies, rejecting algebraic loops.
List<RealAssignment> algebraic_sort(List<RealAssignment> const&);

//! Shared symbolic constructors for continuous, discrete and hybrid systems.
EffectiveVectorMultivariateFunction make_auxiliary_function(Space<Real> const&, List<RealAssignment> const&);
EffectiveVectorMultivariateFunction make_dynamic_function(Space<Real> const&, List<RealAssignment> const&, List<DottedRealAssignment> const&);
EffectiveVectorMultivariateFunction make_reset_function(Space<Real> const&, List<RealAssignment> const&, List<PrimedRealAssignment> const&);
EffectiveScalarMultivariateFunction make_constraint_function(Space<Real> const&, List<RealAssignment> const&, ContinuousPredicate const&, Sign);

} // namespace Ariadne

#endif // ARIADNE_SYSTEM_SYSTEM_FUNCTION_HPP
