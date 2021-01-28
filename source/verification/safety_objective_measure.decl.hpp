/***************************************************************************
 *            verification/safety_objective_measure.decl.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file verification/safety_objective_measure.decl.hpp
 *  \brief Declarations for specific safety objectives
 */

#ifndef ARIADNE_SAFETY_OBJECTIVE_MEASURE_DECL_HPP
#define ARIADNE_SAFETY_OBJECTIVE_MEASURE_DECL_HPP

#include "numeric/dyadic.hpp"
#include "numeric/float_value.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float.decl.hpp"
#include "safety_objective_measure.hpp"

namespace Ariadne {

template class SafetyObjectiveMeasure<PositiveFloatDPUpperBound,Dyadic>;
using TimedRadiusObjective = SafetyObjectiveMeasure<PositiveFloatDPUpperBound,Dyadic>;

} // namespace Ariadne

#endif // ARIADNE_SAFETY_OBJECTIVE_MEASURE_DECL_HPP
