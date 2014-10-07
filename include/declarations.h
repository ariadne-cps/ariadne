/***************************************************************************
 *            declarations.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file declarations.h
 *  \brief Forward declarations of types and classes.
 */

#ifndef ARIADNE_DECLARATIONS_H
#define ARIADNE_DECLARATIONS_H

#include <iosfwd>
#include <boost/logic/tribool.hpp>

namespace Ariadne {

//! Internal name for output stream.
typedef std::ostream OutputStream;

//! Internal name for void type.
typedef void Void;

//! Internal name for builtin boolean type.
typedef bool Bool;
//! Internal name for three-valued logical type.
typedef boost::tribool Tribool;
//! Internal name for string objects.
typedef std::string String;
//! Internal name for builtin unsigned integers.
typedef unsigned int Nat;
//! Internal name for builtin integers.
typedef int Int;



// Numeric declarations
class Integer;
class Rational;
class Real;

class Float;
typedef Float RawFloat;

class ExactFloat;
class ValidatedFloat;
class UpperFloat;
class LowerFloat;
class ApproximateFloat;
typedef UpperFloat PositiveUpperFloat;
typedef PositiveUpperFloat ErrorFloat;

typedef ApproximateFloat ApproximateNumber;
typedef LowerFloat LowerNumber;
typedef UpperFloat UpperNumber;
typedef ValidatedFloat ValidatedNumber;
typedef Real EffectiveNumber;
typedef ExactFloat ExactNumber;

typedef PositiveUpperFloat PositiveUpperNumber;

typedef ErrorFloat NormType;
typedef ErrorFloat ErrorType;
typedef ExactFloat CoefficientType;

typedef ApproximateNumber ApproximateNormType;
typedef ApproximateNumber ApproximateErrorType;
typedef ApproximateNumber ApproximateCoefficientType;

// Information level declarations
struct ExactTag { };
typedef EffectiveNumber EffectiveTag;
typedef ValidatedNumber ValidatedTag;
typedef ApproximateNumber ApproximateTag;

template<class I> struct CanonicalNumberTypedef;
template<> struct CanonicalNumberTypedef<ExactTag> { typedef ExactNumber Type; };
template<> struct CanonicalNumberTypedef<EffectiveTag> { typedef EffectiveNumber Type; };
template<> struct CanonicalNumberTypedef<ValidatedTag> { typedef ValidatedNumber Type; };
template<> struct CanonicalNumberTypedef<ApproximateTag> { typedef ApproximateNumber Type; };
template<class I> using CanonicalNumber = typename CanonicalNumberTypedef<I>::Type;

// Concrete class declarations
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class Vector<Differential<X>>;
template<class X> class Series;

template<class X> class AffineModel;
template<class X> class TaylorModel;
template<class X> class Formula;
template<class X> class Algebra;

template<class X> class Point;
class ExactInterval;
class UpperInterval;
class ApproximateInterval;
class ExactBox;
class UpperBox;
class ApproximateBox;

typedef Vector<RawFloat> RawFloatVector;
typedef Vector<ApproximateFloat> ApproximateFloatVector;
typedef Vector<ValidatedFloat> ValidatedFloatVector;
typedef Vector<ExactFloat> ExactFloatVector;

typedef Vector<ExactInterval> ExactIntervalVector;
typedef Vector<UpperInterval> UpperIntervalVector;
typedef Vector<ApproximateInterval> ApproximateIntervalVector;

typedef Vector<ApproximateNumber> ApproximateVector;
typedef Vector<ValidatedNumber> ValidatedVector;
typedef Vector<EffectiveNumber> EffectiveVector;
typedef Vector<ExactNumber> ExactVector;

typedef Point<ApproximateNumber> ApproximatePoint;
typedef Point<ValidatedNumber> ValidatedPoint;
typedef Point<EffectiveNumber> EffectivePoint;
typedef Point<ExactNumber> ExactPoint;

// Function interface declarations
template<class X> class ScalarFunctionInterface;
template<class X> class VectorFunctionInterface;

typedef ScalarFunctionInterface<ApproximateTag> ApproximateScalarFunctionInterface;
typedef ScalarFunctionInterface<ValidatedTag> ValidatedScalarFunctionInterface;
typedef ScalarFunctionInterface<EffectiveTag> EffectiveScalarFunctionInterface;

typedef VectorFunctionInterface<ApproximateTag> ApproximateVectorFunctionInterface;
typedef VectorFunctionInterface<ValidatedTag> ValidatedVectorFunctionInterface;
typedef VectorFunctionInterface<EffectiveTag> EffectiveVectorFunctionInterface;

// Function declarations
template<class X> class ScalarFunction;
template<class X> class VectorFunction;

typedef ScalarFunction<ApproximateTag> ApproximateScalarFunction;
typedef ScalarFunction<ValidatedTag> ValidatedScalarFunction;
typedef ScalarFunction<EffectiveTag> EffectiveScalarFunction;
typedef EffectiveScalarFunction RealScalarFunction;

typedef VectorFunction<ApproximateTag> ApproximateVectorFunction;
typedef VectorFunction<ValidatedTag> ValidatedVectorFunction;
typedef VectorFunction<EffectiveTag> EffectiveVectorFunction;
typedef EffectiveVectorFunction RealVectorFunction;

// Function model declarations
template<class X> class ScalarFunctionModel;
template<class X> class VectorFunctionModel;

typedef ScalarFunctionModel<ApproximateTag> ApproximateScalarFunctionModel;
typedef ScalarFunctionModel<ValidatedTag> ValidatedScalarFunctionModel;

typedef VectorFunctionModel<ApproximateTag> ApproximateVectorFunctionModel;
typedef VectorFunctionModel<ValidatedTag> ValidatedVectorFunctionModel;


} // namespace Ariadne

#endif
