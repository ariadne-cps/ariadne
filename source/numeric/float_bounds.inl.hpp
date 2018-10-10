/***************************************************************************
 *            float_bounds.inl.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

namespace Ariadne {

template<class F> Bounds<F> Operations<Bounds<F>>::_mul(Bounds<F> const& x1, Bounds<F> const& x2) {
    const F& x1l=x1._l; const F& x1u=x1._u;
    const F& x2l=x2._l; const F& x2u=x2._u;
    F rl,ru;
    if(x1l>=0) {
        if(x2l>=0) {
            rl=mul(down,x1l,x2l); ru=mul(up,x1u,x2u);
        } else if(x2u<=0) {
            rl=mul(down,x1u,x2l); ru=mul(up,x1l,x2u);
        } else {
            rl=mul(down,x1u,x2l); ru=mul(up,x1u,x2u);
        }
    }
    else if(x1u<=0) {
        if(x2l>=0) {
            rl=mul(down,x1l,x2u); ru=mul(up,x1u,x2l);
        } else if(x2u<=0) {
            rl=mul(down,x1u,x2u); ru=mul(up,x1l,x2l);
        } else {
            rl=mul(down,x1l,x2u); ru=mul(up,x1l,x2l);
        }
    } else {
        if(x2l>=0) {
            rl=mul(down,x1l,x2u); ru=mul(up,x1u,x2u);
        } else if(x2u<=0) {
            rl=mul(down,x1u,x2l); ru=mul(up,x1l,x2l);
        } else {
            rl=min(mul(down,x1u,x2l),mul(down,x1l,x2u));
            ru=max(mul(up,x1l,x2l),mul(up,x1u,x2u));
        }
    }
    return Bounds<F>(rl,ru);
}

template<class F> Bounds<F> Operations<Bounds<F>>::_div(Bounds<F> const& x1, Bounds<F> const& x2) {
    const F& x1l=x1.lower_raw(); const F& x1u=x1.upper_raw();
    const F& x2l=x2.lower_raw(); const F& x2u=x2.upper_raw();
    F rl,ru;

    // IMPORTANT: Need to be careful when one of the bounds is 0, since if x2l=-0.0 and x1u>0, then x2l>=0 but x1u/x2l=-inf
    if(x2l>0) {
        if(x1l>=0) {
            rl=div(down,x1l,x2u); ru=div(up,x1u,x2l);
        } else if(x1u<=0) {
            rl=div(down,x1l,x2l); ru=div(up,x1u,x2u);
        } else {
            rl=div(down,x1l,x2l); ru=div(up,x1u,x2l);
        }
    }
    else if(x2u<0) {
        if(x1l>=0) {
            rl=div(down,x1u,x2u); ru=div(up,x1l,x2l);
        } else if(x1u<=0) {
            rl=div(down,x1u,x2l); ru=div(up,x1l,x2u);
        } else {
            rl=div(down,x1u,x2u); ru=div(up,x1l,x2u);
        }
    }
    else {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
        PR pr=max(x1.precision(),x2.precision());
        rl=-F::inf(pr);
        ru=+F::inf(pr);
    }
    return Bounds<F>(rl,ru);
}

template<class F> Bounds<F> Operations<Bounds<F>>::_pi(PR pr) {
    return Bounds<F>(F::pi(down,pr),F::pi(up,pr));
}

template<class F> Bounds<F> Operations<Bounds<F>>::_sin(Bounds<F> const& x) {
    return cos(x-hlf(_pi(x.precision())));
}

template<class F> Bounds<F> Operations<Bounds<F>>::_cos(Bounds<F> const& x) {
    ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
    typename F::RoundingModeType rnd = F::get_rounding_mode();
    PR prec=x.precision();

    const F one(1,prec);
    const Value<F> two(2,prec);
    const Bounds<F> pi_val=_pi(prec);
    const Bounds<F> two_pi_val=2*pi_val;
    if(x.error().raw()>two_pi_val.lower().raw()) { return Bounds<F>(-one,+one); }

    Value<F> n(round(div(near,x.value_raw(),(two_pi_val.value_raw()))));
    Bounds<F> y=x-two*(n*pi_val);

    ARIADNE_ASSERT(y.lower_raw()<=pi_val.upper_raw());
    ARIADNE_ASSERT(y.upper_raw()>=-pi_val.upper_raw());

    F rl,ru;
    if(y.lower_raw()<=-pi_val.lower_raw()) {
        if(y.upper_raw()<=0.0) { rl=-one; ru=cos(up,y.upper_raw()); }
        else { rl=-one; ru=+one; }
    } else if(y.lower_raw()<=0.0) {
        if(y.upper_raw()<=0.0) { rl=cos(down,y.lower_raw()); ru=cos(up,y.upper_raw()); }
        else if(y.upper_raw()<=pi_val.lower_raw()) { rl=cos(down,max(-y.lower_raw(),y.upper_raw())); ru=+one; }
        else { rl=-one; ru=+one; }
    } else if(y.lower_raw()<=pi_val.upper_raw()) {
        if(y.upper_raw()<=pi_val.lower_raw()) { rl=cos(down,y.upper_raw()); ru=cos(up,y.lower_raw()); }
        else if(y.upper_raw()<=two_pi_val.lower_raw()) { rl=-one; ru=cos(up,min(y.lower_raw(),sub(down,two_pi_val.lower_raw(),y.upper_raw()))); }
        else { rl=-one; ru=+one; }
    } else {
        assert(false);
    }

    F::set_rounding_mode(rnd);
    return Bounds<F>(rl,ru);
}


} // namespace Ariadne
