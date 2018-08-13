/***************************************************************************
 *            declarations.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file declarations.hpp
 *  \brief Forward declarations of types and classes.
 */

#ifndef ARIADNE_DECLARATIONS_HPP
#define ARIADNE_DECLARATIONS_HPP

#include <iosfwd>

#include "../utility/metaprogramming.hpp"
#include "../utility/typedefs.hpp"

#include "../numeric/paradigm.hpp"
#include "../numeric/logical.decl.hpp"
#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"

#include "../function/function.decl.hpp"

#include "../geometry/interval.decl.hpp"
#include "../geometry/box.decl.hpp"

namespace Ariadne {

typedef unsigned int uint;

//! Internal name for output stream.
typedef OutputStream OutputStream;

//! Internal name for void type.
typedef void Void;

//! Internal name for builtin boolean type.
typedef bool Bool;
//! Internal name for builtin integers.
typedef int Int;
//! Internal name for builtin unsigned integers.
typedef uint Nat;

// A class containing an exact double-precision value
class ExactDouble;

// Define as a class for consistency with other value types
class String;


typedef SizeType DimensionType;

typedef FloatDPError ValidatedNormType; // FIXME: Remove this typedef
typedef FloatDPApproximation ApproximateNormType; // FIXME: Remove this typedef

typedef FloatDPError NormType; // FIXME: Remove this typedef
typedef FloatDPError ErrorType; // FIXME: Remove this typedef
typedef FloatDPApproximation ApproximateErrorType; // FIXME: Remove this typedef

template<class X> struct InformationTypedef;
template<> struct InformationTypedef<Real> { typedef EffectiveTag Type; };
template<class P> struct InformationTypedef<Number<P>> { typedef P Type; };
template<class X> using InformationTag = typename InformationTypedef<X>::Type;

template<class X> using Scalar = X;
// Concrete class declarations
template<class X> class Vector;
template<class X> class Covector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class UnivariateDifferential;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class Series;

template<class P, class F> class AffineModel;
template<class P, class F> class TaylorModel;
template<class X> class Formula;
template<class X> class Algebra;

template<class X> class Point;

typedef Vector<Rational> RationalVector;
typedef Vector<Real> RealVector;
typedef Vector<FloatDP> FloatVector;
typedef Vector<RawFloatDP> RawFloatVector;
typedef Vector<FloatDPApproximation> FloatApproximationVector;
typedef Vector<FloatDPBounds> FloatBoundsVector;
typedef Vector<FloatDPValue> ExactFloatVector;

typedef Vector<ApproximateNumericType> ApproximateVector;
typedef Vector<ValidatedNumericType> ValidatedVector;
typedef Vector<EffectiveNumericType> EffectiveVector;
typedef Vector<ExactNumericType> ExactVector;

typedef Matrix<Rational> RationalMatrix;
typedef Matrix<Real> RealMatrix;
typedef Matrix<RawFloatDP> RawFloatMatrix;
typedef Matrix<FloatDP> FloatMatrix;
typedef Matrix<FloatDPApproximation> FloatApproximationMatrix;
typedef Matrix<FloatDPBounds> FloatBoundsMatrix;
typedef Matrix<FloatDPValue> ExactFloatMatrix;

typedef Point<ApproximateNumericType> ApproximatePoint;
typedef Point<ValidatedNumericType> ValidatedPoint;
typedef Point<EffectiveNumericType> EffectivePoint;
typedef Point<ExactNumericType> ExactPoint;


} // namespace Ariadne

#endif
