/***************************************************************************
 *            verification/safety_objective_measure.hpp
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

/*! \file verification/safety_objective_measure.hpp
 *  \brief A class for holding the measure of a safety objective
 */

#ifndef ARIADNE_SAFETY_OBJECTIVE_MEASURE_HPP
#define ARIADNE_SAFETY_OBJECTIVE_MEASURE_HPP

#include "symbolic/variable.hpp"
#include "symbolic/expression.decl.hpp"

namespace Ariadne {

//! \brief The safety objective measure requires a measure type M and a time type T
template<class M, class T>
class SafetyObjectiveMeasure {
  public:
    SafetyObjectiveMeasure(RealVariable const& variable, M const& measure, T const& time) : _variable(variable), _measure(measure), _time(time) { }

    RealVariable const& variable() const { return _variable; }
    M const& measure() const { return _measure; }
    T const& time() const { return _time; }
    bool operator<(SafetyObjectiveMeasure<M,T> const& measure) const { return _time < measure._time; }
  private:
    RealVariable const _variable;
    M const _measure;
    T const _time;
};

} // namespace Ariadne

#endif // ARIADNE_SAFETY_OBJECTIVE_MEASURE_HPP
