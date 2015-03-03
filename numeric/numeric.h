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
 *  \brief Number classes. File suitable for use as a pre-compiled header.
 */

#ifndef ARIADNE_NUMERIC_H
#define ARIADNE_NUMERIC_H

#include "utility/standard.h"

#include "config.h"

#include "utility/declarations.h"

#include "numeric/logical.h"
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/decimal.h"
#include "numeric/dyadic.h"
#include "numeric/float.h"
#include "numeric/real.h"
#include "numeric/number.h"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Cast one %Ariadne numerical type or builtin numerical type to another.
template<class R, class A> inline R numeric_cast(const A& a) { return R(a); }
template<> inline Int numeric_cast(const Float64& a) { return Int(a.get_d()); }
template<> inline Int numeric_cast(const FloatMP& a) { return Int(a.get_d()); }
template<> inline double numeric_cast(const Float64& a) { return a.get_d(); }
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline double numeric_cast(const ExactFloat64& a) { return a.get_d(); }
template<> inline double numeric_cast(const BoundedFloat64& a) { return a.get_d(); }
template<> inline double numeric_cast(const ApproximateFloat64& a) { return a.get_d(); }
template<> inline float numeric_cast(const double& a) { return a; }
template<> inline float numeric_cast(const Float64& a) { return a.get_d(); }
template<> inline float numeric_cast(const Real& a) { return a.get_d(); }
template<> inline Float64 numeric_cast(const ExactFloat64& a) { return a.raw(); }

template<> inline Real numeric_cast(const Float64& a) { return Real(ExactFloat64(a)); }
template<> inline Real numeric_cast(const ExactFloat64& a) { return Real(a); }
template<> inline Real numeric_cast(const BoundedFloat64& a) { return Real(cast_exact(ApproximateFloat64(a))); }

} // namespace Ariadne

#endif
