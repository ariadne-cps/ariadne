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

#include "numeric/paradigm.h"
#include "is_number.h"

namespace Ariadne {

/************ Number *********************************************************/

class DivideByZeroError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

template<class X> struct IsNumber;
template<class P=Void> class Number;

class Nat32;
class Nat64;
class Int32;
class Int64;

class Dyadic;
class Decimal;

class Integer;
class Rational;
class Real;

template<> struct IsNumber<Integer>;
template<> struct IsNumber<Rational>;
template<> struct IsNumber<Real>;

/*
using ExactNumber=Number<Exact>;
using EffectiveNumber=Number<Effective>;
using ValidatedNumber=Number<Validated>;
using UpperNumber=Number<Upper>;
using LowerNumber=Number<Lower>;
using ApproximateNumber=Number<Approximate>;
*/
}

#include "float.decl.h"

namespace Ariadne {

using ExactNumber=ExactFloat;
using EffectiveNumber=Real;
using ValidatedNumber=ValidatedFloat;
using UpperNumber=UpperFloat;
using LowerNumber=LowerFloat;
using ApproximateNumber=ApproximateFloat;
using PositiveUpperNumber=PositiveUpperFloat;

} // namespace Ariadne

#endif
