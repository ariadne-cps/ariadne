/***************************************************************************
 *            numeric/float_factory.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file numeric/float_factory.hpp
 *  \brief Factories for creating floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_FACTORY_HPP
#define ARIADNE_FLOAT_FACTORY_HPP

#include "float.decl.hpp"

namespace Ariadne {

template<class PR> class FloatFactory {
    PR _pr; typedef RawFloatType<PR> F;
  public:
    typedef PR PrecisionType;
    typedef PR PropertiesType;
    FloatFactory(PR const& pr) : _pr(pr) { }
    PR precision() const { return this->_pr; }
    PR properties() const { return this->_pr; }
  public:
    FloatApproximation<PR> create(Number<ApproximateTag> const& y);
    FloatLowerBound<PR> create(Number<ValidatedLowerTag> const& y);
    FloatUpperBound<PR> create(Number<ValidatedUpperTag> const& y);
    FloatBounds<PR> create(Number<ValidatedTag> const& y);
    FloatBounds<PR> create(Number<EffectiveTag> const& y);
    FloatBounds<PR> create(Number<ExactTag> const& y);
    FloatBounds<PR> create(Real const& y);
    FloatBounds<PR> create(Rational const& y);
    FloatBounds<PR> create(Dyadic const& y);
    FloatBounds<PR> create(Integer const& y);
    FloatValue<PR> create(ExactDouble const& y);
    FloatValue<PR> create(Dyadic const& y, ExactTag);
    FloatValue<PR> create(Integer const& y, ExactTag);
    template<class N, EnableIf<IsBuiltinSignedIntegral<N>> =dummy> FloatValue<PR> create(N const& y);
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR> create(M const& y);
    template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> FloatApproximation<PR> create(D const& y);
    PositiveFloatApproximation<PR> create(PositiveNumber<ApproximateTag> const& y);
    PositiveFloatLowerBound<PR> create(PositiveNumber<ValidatedLowerTag> const& y);
    PositiveFloatUpperBound<PR> create(PositiveNumber<ValidatedUpperTag> const& y);
    PositiveFloatBounds<PR> create(PositiveNumber<ValidatedTag> const& y);
};

template<class Y, class PR> using ConcreteType = decltype(declval<FloatFactory<PR>>().create(declval<Y>()));
template<class Y, class PR> inline decltype(auto) make_float(Y const& y, PR pr) { return float_factory(pr).create(y); }

template<class PR> inline FloatFactory<PR> float_factory(PR pr) { return FloatFactory<PR>(pr); }
template<class F> inline FloatFactory<PrecisionType<F>> factory(Approximation<F> const& flt);
template<class F> inline FloatFactory<PrecisionType<F>> factory(LowerBound<F> const& flt);
template<class F> inline FloatFactory<PrecisionType<F>> factory(UpperBound<F> const& flt);
template<class F> inline FloatFactory<PrecisionType<F>> factory(Bounds<F> const& flt);
template<class F, class FE> inline FloatFactory<PrecisionType<F>> factory(Ball<F,FE> const& flt);
template<class F> inline FloatFactory<PrecisionType<F>> factory(Value<F> const& flt);

/*
template<class F> inline FloatFactory<PR> factory(Approximation<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(LowerBound<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(UpperBound<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(Bounds<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F, class FE> inline FloatFactory<PR> factory(FloatBall<PR,PRE> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(Value<F> const& flt) { return FloatFactory<PR>(flt.precision()); }

template<class F> inline Approximation<F> FloatFactory<PR>::(Number<ApproximateTag> const& y) { return Approximation<F>(y,_pr); }
template<class F> inline LowerBound<F> FloatFactory<PR>::(Number<ValidatedLowerTag> const& y) { return LowerBound<F>(y,_pr); }
template<class F> inline UpperBound<F> FloatFactory<PR>::(Number<ValidatedUpperTag> const& y) { return UpperBound<F>(y,_pr); }

template<class F> inline Bounds<F> FloatFactory<PR>::(Number<ValidatedTag> const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Number<EffectiveTag> const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Number<ExactTag> const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Real const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Rational const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Dyadic const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Integer const& y) { return Bounds<F>(y,_pr); }

template<class F> inline Value<F> FloatFactory<PR>::(Dyadic const& y, ExactTag) { return Value<F>(y,_pr); }
template<class F> inline Value<F> FloatFactory<PR>::(Integer const& y, ExactTag) { return Value<F>(y,_pr); }

template<class F> inline template<class N, EnableIf<IsBuiltinSignedIntegral<N>> =dummy> Value<F> FloatFactory<PR>::(N const& y) { return Value<F>(y,_pr); }
template<class F> inline template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> PositiveValue<F> FloatFactory<PR>::(M const& y) { return PositiveValue<F>(y,_pr); }
template<class F> inline template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> Approximation<F> FloatFactory<PR>::(D const& y) { return Approximation<F>(F(y,_pr)); }
*/

} // namespace Ariadne

#endif
