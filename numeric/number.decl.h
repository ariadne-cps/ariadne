/***************************************************************************
 *            numeric/number.decl.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/number.decl.h
 *  \brief
 */

#ifndef ARIADNE_NUMBER_DECL_H
#define ARIADNE_NUMBER_DECL_H

#include <stdexcept>

#include "utility/metaprogramming.h"
#include "utility/typedefs.h"
#include "numeric/paradigm.h"

namespace Ariadne {

/************ Number *********************************************************/

class DivideByZeroError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

template<class X> struct IsNumericType : False { };
template<class X, class T> using EnableIfNumericType = EnableIf<IsNumericType<X>,T>;

typedef uint Nat;
typedef int Int;
typedef double Dbl;

class Nat32;
class Nat64;
class Int32;
class Int64;

class Dyadic;
class Decimal;

class Integer;
class Rational;
class Real;

template<class P=Void> class Number;

template<> struct IsNumericType<Nat>;
template<> struct IsNumericType<Int>;
template<> struct IsNumericType<Dbl>;

template<> struct IsNumericType<Integer>;
template<> struct IsNumericType<Rational>;
template<> struct IsNumericType<Real>;

template<class X> struct IsGenericNumericType : IsConvertible<X,Real> { };
template<> struct IsGenericNumericType<Real> : True { };
template<> struct IsGenericNumericType<Dbl> : True { };
template<class P> struct IsGenericNumericType<Number<P>> : True { };


using ExactNumber=Number<ExactTag>;
using EffectiveNumber=Number<EffectiveTag>;
using ValidatedNumber=Number<ValidatedTag>;
using UpperNumber=Number<UpperTag>;
using LowerNumber=Number<LowerTag>;
using ApproximateNumber=Number<ApproximateTag>;
using PositiveUpperNumber=Number<PositiveUpperTag>;

} // namespace Ariadne

#include "float.decl.h"

namespace Ariadne {

using ExactNumericType=ExactFloat64;
using EffectiveNumericType=Real;
using ValidatedNumericType=BoundedFloat64;
using UpperNumericType=UpperFloat64;
using LowerNumericType=LowerFloat64;
using ApproximateNumericType=ApproximateFloat64;
using PositiveUpperNumericType=PositiveUpperFloat64;

} // namespace Ariadne

#endif
