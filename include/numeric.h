/***************************************************************************
 *            numeric.h
 *
 *  Copyright 2008-10  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file numeric.h
 *  \brief Numerical classes.
 */
#ifndef ARIADNE_NUMERIC_H
#define ARIADNE_NUMERIC_H

#include "integer.h"
#include "rational.h"
#include "float.h"
#include "interval.h"
#include "real.h"

namespace Ariadne {

template<class R, class A> inline R numeric_cast(const A& a);
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline Real numeric_cast(const Float& a) { return Real(a); }
template<> inline Real numeric_cast(const Interval& a) { return Real(a); }

//! \ingroup NumericModule \related Float \related Interval \related Real
//! \brief Cast one %Ariadne numerical type or builtin numerical type to another.
template<class R, class A> inline R numeric_cast(const A& a) { return R(a); }
template<> inline int numeric_cast(const Float& a) { return int(a.get_d()); }
template<> inline double numeric_cast(const Float& a) { return a.get_d(); }
template<> inline double numeric_cast(const Interval& a) { return a.get_d(); }
template<> inline Float numeric_cast(const Interval& a) { return a.midpoint(); }
template<> inline Interval numeric_cast(const Float& a) { return Interval(a); }

//! \ingroup NumericModule \related Float
//! \brief Converts \a e to an object of type \a X, which may either be an
//! \c Float or \c Interval, with the semantics that \a e denotes and error bound.
//! Returns the float 0.0 (since floating-point computations do not keep track of errors)
//! and the interval [-e,+e].
template<class X> inline X convert_error(const Float& e);
template<> inline Float convert_error<Float>(const Float& e) { return 0.0; }
template<> inline Interval convert_error<Interval>(const Float& e) { return Interval(-e,+e); }

// Use 'enable_if' style template to restrict allowable instances. See the Boost documentation
// for enable_if to see how this works.
template<class X, class T> struct enable_if_numeric { };
template<class T> struct enable_if_numeric<unsigned int,T> { typedef T type; };
template<class T> struct enable_if_numeric<int,T> { typedef T type; };
template<class T> struct enable_if_numeric<double,T> { typedef T type; };
template<class T> struct enable_if_numeric<Float,T> { typedef T type; };
template<class T> struct enable_if_numeric<Interval,T> { typedef T type; };
template<class T> struct enable_if_numeric<Real,T> { typedef T type; };


} // namespace Ariadne

#endif
