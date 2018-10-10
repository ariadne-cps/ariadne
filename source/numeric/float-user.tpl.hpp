/***************************************************************************
 *            float-user.tpl.hpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "../config.hpp"
#include "../utility/macros.hpp"
#include "../utility/exceptions.hpp"

#include "../numeric/float-user.hpp"

#include "../numeric/integer.hpp"
#include "../numeric/dyadic.hpp"
#include "../numeric/decimal.hpp"
#include "../numeric/rational.hpp"
#include "../numeric/real.hpp"

#include "../numeric/number_wrapper.hpp"

namespace Ariadne {

template<class PR> Approximation<RawFloatType<PR>> make_float(Number<ApproximateTag> x, PR pr) { return Approximation<RawFloatType<PR>>(x,pr); }
template<class PR> LowerBound<RawFloatType<PR>> make_float(Number<ValidatedLowerTag> x, PR pr) { return LowerBound<RawFloatType<PR>>(x,pr); }
template<class PR> UpperBound<RawFloatType<PR>> make_float(Number<ValidatedUpperTag> x, PR pr) { return UpperBound<RawFloatType<PR>>(x,pr); }
template<class PR> Bounds<RawFloatType<PR>> make_float(Number<ValidatedTag> x, PR pr) { return Bounds<RawFloatType<PR>>(x,pr); }
template<class PR, class PRE> Ball<RawFloatType<PR>,RawFloatType<PRE>> make_float(Number<ValidatedTag> x, PR pr, PRE pre) { return Ball<RawFloatType<PR>,RawFloatType<PRE>>(x,pr,pre); }
template<class PR> Bounds<RawFloatType<PR>> make_float(Number<EffectiveTag> x, PR pr) { return Bounds<RawFloatType<PR>>(x,pr); }
template<class PR> Bounds<RawFloatType<PR>> make_float(Number<ExactTag> x, PR pr) { return Bounds<RawFloatType<PR>>(x,pr); }
template<class PR> Bounds<RawFloatType<PR>> make_float(Real r) { return Bounds<RawFloatType<PR>>(r); }
template<class PR> Bounds<RawFloatType<PR>> make_float(Rational q) { return Bounds<RawFloatType<PR>>(q); }
template<class PR> Value<RawFloatType<PR>> make_float(Integer z) { return Value<RawFloatType<PR>>(z); }


template<class F> Approximation<F>::Approximation(LowerBound<F> const& x) : _a(x.raw()) { }
template<class F> Approximation<F>::Approximation(UpperBound<F> const& x) : _a(x.raw()) { }
template<class F> Approximation<F>::Approximation(Bounds<F> const& x) : _a(x.value_raw()) { }
template<class F> Approximation<F>::Approximation(Value<F> const& x) : _a(x.raw()) { }
template<class F> Approximation<F>::Approximation(Error<F> const& x) : _a(x.raw()) { }

template<class F> LowerBound<F>::LowerBound(Bounds<F> const& x) : _l(x.lower_raw()) { }
template<class F> LowerBound<F>::LowerBound(Ball<F> const& x) : _l(x.lower_raw()) { }
template<class F> LowerBound<F>::LowerBound(Value<F> const& x) : _l(x.raw()) { }

template<class F> UpperBound<F>::UpperBound(Bounds<F> const& x) : _u(x.upper_raw()) { }
template<class F> UpperBound<F>::UpperBound(Ball<F> const& x) : _u(x.upper_raw()) { }
template<class F> UpperBound<F>::UpperBound(Value<F> const& x) : _u(x.raw()) { }
template<class F> UpperBound<F>::UpperBound(Error<F> const& x) : _u(x.raw()) { }

template<class F> Bounds<F>::Bounds(Value<F> const& x) : _l(x.raw()), _u(x.raw()) { }

template<class F, class FE> Ball<F,FE>::Ball(Bounds<F> const& x)
    : _v(x.value_raw()) , _e(_make_error<FE>(x.error_raw())) { }

template<class F, class FE> Ball<F,FE>::Ball(Value<F> const& x)
    : _v(x.raw()), _e(_error_precision<PRE>(x.precision())) { }


template<class F, class FE> Ball<F,FE> Ball<F,FE>::pm(Error<FE> const& e) const {
    return Ball<F,FE>(this->_v,add(up,this->_e,e._e)); }

template<class F> Bounds<F> Bounds<F>::pm(Error<F> e) const {
    return Bounds<F>(sub(down,this->_l,e._e),add(up,this->_u,e._e)); }


} // namespace Ariadne
