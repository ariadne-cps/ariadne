/***************************************************************************
 *            dynamics/vector_field.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "function/functional.hpp"
#include "config.hpp"

#include "utility/macros.hpp"
#include "utility/array.hpp"
#include "utility/tuple.hpp"
#include "utility/stlio.hpp"
#include "utility/container.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"

#include "solvers/integrator.hpp"

#include "output/logging.hpp"

#include "dynamics/vector_field.hpp"

#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"

namespace Ariadne {

EffectiveVectorMultivariateFunction make_auxiliary_function(
    Space<Real> const& state_space,
    List<RealAssignment> const& algebraic);

EffectiveVectorMultivariateFunction make_dynamic_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<DottedRealAssignment> const& differential);


VectorField::VectorField(List<DottedRealAssignment> const& dynamics)
    : VectorField(dynamics, List<RealAssignment>())
{
}

VectorField::VectorField(List<DottedRealAssignment> const& dynamics, List<RealAssignment> const& auxiliary)
    : _dynamics(dynamics), _auxiliary(auxiliary)
    , _dynamic_function(make_dynamic_function(left_hand_sides(dynamics),auxiliary,dynamics))
    , _auxiliary_function(make_auxiliary_function(left_hand_sides(dynamics),auxiliary))
{
}

VectorField::VectorField(EffectiveVectorMultivariateFunction const& function)
    : _dynamic_function(function), _auxiliary_function(0u,function.domain())
{
    ARIADNE_PRECONDITION(function.result_size()==function.argument_size());
}

RealSpace VectorField::state_space() const {
    return RealSpace(left_hand_sides(this->_dynamics));
}

RealSpace VectorField::auxiliary_space() const {
    return RealSpace(left_hand_sides(this->_auxiliary));
}


OutputStream& operator<<(OutputStream& os, const VectorField& vf) {
    os << "VectorField( dynamic_function = " << vf.dynamic_function() << ", "
          "auxiliary_function = " << vf.auxiliary_function() << ", "
          "dynamics = " << vf._dynamics << ", "
          "auxiliary = " << vf._auxiliary << ")";
    return os;
}


}  // namespace Ariadne

