/***************************************************************************
 *            numeric/number.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/number.cpp
 *  \brief
 */



#include "../utility/module.hpp"
#include "../numeric/paradigm.hpp"

#include "casts.hpp"
#include "number.hpp"
#include "logical.hpp"
#include "integer.hpp"
#include "dyadic.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"

#include "number_interface.hpp"
#include "number_wrapper.hpp"

/************ Number *********************************************************/

namespace Ariadne {

using NumberHandle = Handle<NumberInterface>;

// Declare operations for Integer and Rational numbers
Real sqrt(Real const&);
Real exp(Real const&);
Real log(Real const&);
Real sin(Real const&);
Real cos(Real const&);
Real tan(Real const&);
Real asin(Real const&);
Real acos(Real const&);
Real atan(Real const&);

template<> struct DispatchingTraits<Integer> { typedef Aware<Integer> AwareOfTypes; };
template class NumberWrapper<Integer>;
template<> struct DispatchingTraits<Dyadic> { typedef Aware<Dyadic,Integer> AwareOfTypes; };
template class NumberWrapper<Dyadic>;
template<> struct DispatchingTraits<Rational> { typedef Aware<Rational,Dyadic,Integer> AwareOfTypes; };
template class NumberWrapper<Rational>;

template<> struct DispatchingTraits<Real> { typedef Aware<Real> AwareOfTypes; };
template class NumberWrapper<Real>;

ExactDouble::operator ExactNumber() const { return Dyadic(*this).operator ExactNumber(); }
//ExactDouble::operator ExactNumber() const { if(std::isfinite(this->get_d()) { return Dyadic(*this).operator ExactNumber(); } else { return ExactNumber(new NumberWrapper<ExactDouble>(*this)); } }
Integer::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Integer>(*this)); }
Dyadic::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Dyadic>(*this)); }
Rational::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Rational>(*this)); }
Real::operator EffectiveNumber() const { return EffectiveNumber(new NumberWrapper<Real>(*this)); }

FloatDPBall NumberInterface::_get(MetricTag p, DoublePrecision pr) const { return this->_get(p,pr,pr); }
FloatMPBall NumberInterface::_get(MetricTag p, MultiplePrecision pr) const { return this->_get(p,pr,pr); }

//template<class X> struct DispatchingTraits<Value<X>> { typedef Aware<Value<X>,Integer,Dyadic> AwareOfTypes; };
template<class F> struct DispatchingTraits<Value<F>> {
    typedef Aware<Dyadic,Integer> AwareOfTypes; };
template<class F> struct DispatchingTraits<Bounds<F>> {
    typedef Aware<Value<F>,Bounds<F>,Integer,Dyadic,Rational> AwareOfTypes; };
template<class F> struct DispatchingTraits<Approximation<F>> {
    typedef Aware<Value<F>,Bounds<F>,Approximation<F>,Integer,Dyadic,Rational> AwareOfTypes; };

template<class F, class FE> struct DispatchingTraits<Ball<F,FE>> {
    typedef Aware<Ball<F,FE>,Value<F>,Bounds<F>,Approximation<F>,Integer,Dyadic,Rational> AwareOfTypes; };

template class NumberWrapper<FloatDPApproximation>;
//template class NumberWrapper<FloatDPLowerBound>;
//template class NumberWrapper<FloatDPUpperBound>;
template class NumberWrapper<FloatDPBounds>;
template class NumberWrapper<FloatDPBall>;
template class NumberWrapper<FloatDPValue>;

DyadicBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatDPBounds>(FloatDPBounds(*this,dp))); }

template<> FloatDPApproximation::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<FloatDPApproximation>(*this)); }
//template<> FloatDPLowerBound::operator ValidatedLowerNumber() const { return ValidatedLowerNumber(new NumberWrapper<FloatDPLowerBound>(*this)); }
//template<> FloatDPUpperBound::operator ValidatedUpperNumber() const { return ValidatedUpperNumber(new NumberWrapper<FloatDPUpperBound>(*this)); }
template<> FloatDPBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatDPBounds>(*this)); }
//template<> FloatDPBall::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatDPBall>(*this)); }
template<> FloatDPValue::operator ExactNumber() const { return ExactNumber(new NumberWrapper<FloatDPValue>(*this)); }

//template<> FloatDPError::operator ValidatedErrorNumber() const { return ValidatedErrorNumber(new NumberWrapper<FloatDPError>(*this)); }
//template<> FloatMPError::operator ValidatedErrorNumber() const { return ValidatedErrorNumber(new NumberWrapper<FloatMPError>(*this)); }
template<> FloatDPError::operator ValidatedErrorNumber() const { return ValidatedErrorNumber(new NumberWrapper<FloatDPValue>(cast_exact(*this))); }
template<> FloatMPError::operator ValidatedErrorNumber() const { return ValidatedErrorNumber(new NumberWrapper<FloatMPValue>(cast_exact(*this))); }

template class NumberWrapper<FloatMPApproximation>;
//template class NumberWrapper<FloatMPLowerBound>;
//template class NumberWrapper<FloatMPUpperBound>;
template class NumberWrapper<FloatMPBounds>;
template class NumberWrapper<FloatMPBall>;
template class NumberWrapper<FloatMPValue>;

template<> FloatMPApproximation::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<FloatMPApproximation>(*this)); }
//template<> FloatMPLowerBound::operator ValidatedLowerNumber() const { return ValidatedLowerNumber(new NumberWrapper<FloatMPLowerBound>(*this)); }
//template<> FloatMPUpperBound::operator ValidatedUpperNumber() const { return ValidatedUpperNumber(new NumberWrapper<FloatMPUpperBound>(*this)); }
template<> FloatMPBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatMPBounds>(*this)); }
//template<> FloatMPBall::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatMPBall>(*this)); }
template<> FloatMPValue::operator ExactNumber() const { return ExactNumber(new NumberWrapper<FloatMPValue>(*this)); }

template<> String class_name<NumberHandle>() { return "NumberHandle"; }

inline Bool refines(Number<UpperTag> const& y1, Number<UpperTag> const& y2) {
    return y1.get(dp).raw() <= y2.get(dp).raw(); }

Positive<ValidatedUpperNumber> mag(Positive<ValidatedUpperNumber> const& y) {
    return y; }


template<> String class_name<ApproximateNumber>() { return "ApproximateNumber"; }
template<> String class_name<ValidatedLowerNumber>() { return "ValidatedLowerNumber"; }
template<> String class_name<ValidatedUpperNumber>() { return "ValidatedUpperNumber"; }
template<> String class_name<ValidatedNumber>() { return "ValidatedNumber"; }
template<> String class_name<EffectiveNumber>() { return "EffectiveNumber"; }
template<> String class_name<ExactNumber>() { return "ExactNumber"; }

} // namespace Ariadne
