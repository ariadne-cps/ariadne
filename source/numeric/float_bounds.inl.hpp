/***************************************************************************
 *            float_bounds.inl.hpp
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


#ifndef ARIADNE_FLOAT_BOUNDS_INL_HPP
#define ARIADNE_FLOAT_BOUNDS_INL_HPP

#include "float_bounds.hpp"
#include "float_value.hpp"

namespace Ariadne {

template<class F> inline auto Operations<Bounds<F>>::_add(Bounds<F> const& x1, Value<F> const& x2) -> Bounds<F> {
    return Bounds<F>(add(down,x1._l,x2._v),add(up,x1._u,x2._v));
}

template<class F> inline auto Operations<Bounds<F>>::_add(Value<F> const& x1, Bounds<F> const& x2) -> Bounds<F> {
    return Bounds<F>(add(down,x1._v,x2._l),add(down,x1._v,x2._u));
}

template<class F> inline auto Operations<Bounds<F>>::_sub(Bounds<F> const& x1, Value<F> const& x2) -> Bounds<F> {
    return Bounds<F>(sub(down,x1._l,x2._v),sub(up,x1._u,x2._v));
}

template<class F> inline auto Operations<Bounds<F>>::_sub(Value<F> const& x1, Bounds<F> const& x2) -> Bounds<F> {
    return Bounds<F>(sub(down,x1._v,x2._u),sub(up,x1._v,x2._l));
}

template<class F> inline auto Operations<Bounds<F>>::_mul(Bounds<F> const& x1, Value<F> const& x2) -> Bounds<F>
{
    const F& x1l=x1.lower_raw(); const F& x1u=x1.upper_raw();
    const F& x2v=x2.raw();
    F rl,ru;
    if(x2v>=0.0) {
        rl=mul(down,x1l,x2v); ru=mul(up,x1u,x2v);
    } else {
        rl=mul(down,x1u,x2v); ru=mul(up,x1l,x2v);
    }
    return Bounds<F>(rl,ru);
}

template<class F> inline auto Operations<Bounds<F>>::_mul(Value<F> const& x1, Bounds<F> const& x2) -> Bounds<F>
{
    const F& x1v=x1.raw();
    const F& x2l=x2.lower_raw(); const F& x2u=x2.upper_raw();
    F rl,ru;
    if(x1v>=0.0) {
        rl=mul(down,x1v,x2l); ru=mul(up,x1v,x2u);
    } else {
        rl=mul(down,x1v,x2u); ru=mul(up,x1v,x2l);
    }
    return Bounds<F>(rl,ru);
}

template<class F> inline auto Operations<Bounds<F>>::_div(Bounds<F> const& x1, Value<F> const& x2) -> Bounds<F>
{
    const F& x1l=x1.lower_raw();
    const F& x1u=x1.upper_raw();
    const F& x2v=x2.raw();
    F rl,ru;
    if(x2v>0) {
        rl=div(down,x1l,x2v); ru=div(up,x1u,x2v);
    } else if(x2v<0) {
        rl=div(down,x1u,x2v); ru=div(up,x1l,x2v);
    } else {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds const& x1, FloatValue x2)","x1="<<x1<<", x2="<<x2);
        auto pr=min(x1.precision(),x2.precision());
        rl=-F::inf(pr);
        ru=+F::inf(pr);
    }
    return Bounds<F>(rl,ru);
}

template<class F> inline auto Operations<Bounds<F>>::_div(Value<F> const& x1, Bounds<F> const& x2) -> Bounds<F>
{
    const F& x1v=x1.raw();
    const F& i2l=x2.lower_raw();
    const F& i2u=x2.upper_raw();
    F rl,ru;
    if(i2l<=0 && i2u>=0) {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatValue const& x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
        auto pr=min(x1.precision(),x2.precision());
        rl=-F::inf(pr);
        ru=+F::inf(pr);
    } else if(x1v>=0) {
        rl=div(down,x1v,i2u); ru=div(up,x1v,i2l);
    } else {
        rl=div(down,x1v,i2l); ru=div(up,x1v,i2u);
    }
    return Bounds<F>(rl,ru);
}



} // namespace Ariadne

#endif
