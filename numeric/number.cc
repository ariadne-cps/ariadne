/***************************************************************************
 *            number.cc
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

/*! \file number.cc
 *  \brief
 */



#include "utility/module.h"
#include "numeric/paradigm.h"

#include "number.h"
#include "logical.h"
#include "integer.h"
#include "rational.h"
#include "real.h"
#include "float64.h"
#include "floatmp.h"

#include "number_interface.h"
#include "number_wrapper.h"

/************ Number *********************************************************/

namespace Ariadne {

using NumberHandle = Handle<NumberInterface>;

// Declare operations for Integer and Rational numbers
Real sqrt(Real);
Real exp(Real);
Real log(Real);
Real sin(Real);
Real cos(Real);
Real tan(Real);
Real atan(Real);

template class NumberWrapper<Integer>;
template class NumberWrapper<Rational>;
template class NumberWrapper<Real>;

Integer::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Integer>(*this)); }
Rational::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Rational>(*this)); }
Real::operator EffectiveNumber() const { return EffectiveNumber(new NumberWrapper<Real>(*this)); }

template class NumberWrapper<ApprxFloat64>;
template class NumberWrapper<LowerFloat64>;
template class NumberWrapper<UpperFloat64>;
template class NumberWrapper<BoundFloat64>;
template class NumberWrapper<MetrcFloat64>;
template class NumberWrapper<ExactFloat64>;

ExactFloat64::operator ExactNumber() const { return ExactNumber(new NumberWrapper<ExactFloat64>(*this)); }
MetrcFloat64::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<MetrcFloat64>(*this)); }
BoundFloat64::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<BoundFloat64>(*this)); }
UpperFloat64::operator UpperNumber() const { return UpperNumber(new NumberWrapper<UpperFloat64>(*this)); }
LowerFloat64::operator LowerNumber() const { return LowerNumber(new NumberWrapper<LowerFloat64>(*this)); }
ApprxFloat64::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<ApprxFloat64>(*this)); }

template class NumberWrapper<ApprxFloatMP>;
//template class NumberWrapper<LowerFloatMP>;
//template class NumberWrapper<UpperFloatMP>;
template class NumberWrapper<BoundFloatMP>;
template class NumberWrapper<MetrcFloatMP>;
//template class NumberWrapper<ExactFloatMP>;

template<> String class_name<NumberHandle>() { return "NumberHandle"; }

} // namespace Ariadne
