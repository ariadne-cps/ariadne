/***************************************************************************
 *            numeric.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file numeric.hpp
 *  \brief Number classes. File suitable for use as a pre-compiled header.
 */

#ifndef ARIADNE_NUMERIC_HPP
#define ARIADNE_NUMERIC_HPP

#include "../config.hpp"

#include "../utility/standard.hpp"
#include "../utility/declarations.hpp"

#include "../numeric/logical.hpp"
#include "../numeric/integer.hpp"
#include "../numeric/rational.hpp"
#include "../numeric/decimal.hpp"
#include "../numeric/dyadic.hpp"
#include "../numeric/float.hpp"
#include "../numeric/real.hpp"
#include "../numeric/number.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Cast one %Ariadne numerical type or builtin numerical type to another.
template<class R, class A> inline R numeric_cast(const A& a) { return R(a); }
template<> inline Int numeric_cast(const FloatDP& a) { return Int(a.get_d()); }
template<> inline Int numeric_cast(const FloatMP& a) { return Int(a.get_d()); }
template<> inline double numeric_cast(const FloatDP& a) { return a.get_d(); }
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline double numeric_cast(const FloatDPValue& a) { return a.get_d(); }
template<> inline double numeric_cast(const FloatDPBounds& a) { return a.get_d(); }
template<> inline double numeric_cast(const FloatDPApproximation& a) { return a.get_d(); }
template<> inline float numeric_cast(const double& a) { return a; }
template<> inline float numeric_cast(const FloatDP& a) { return a.get_d(); }
template<> inline float numeric_cast(const Real& a) { return a.get_d(); }
template<> inline FloatDP numeric_cast(const FloatDPValue& a) { return a.raw(); }

template<> inline Real numeric_cast(const FloatDP& a) { return Real(ExactDouble(a.get_d())); }
template<> inline Real numeric_cast(const FloatDPValue& a) { return numeric_cast<Real>(a.raw()); }
template<> inline Real numeric_cast(const FloatDPBounds& a) { return numeric_cast<Real>(FloatDPApproximation(a).raw()); }

template<> inline FloatDPBall numeric_cast(const Real& a) { return FloatDPBall(a,dp); }
template<> inline FloatDPBounds numeric_cast(const Real& a) { return FloatDPBounds(a,dp); }
template<> inline FloatDPApproximation numeric_cast(const Real& a) { return FloatDPApproximation(a,dp); }

} // namespace Ariadne

#endif
