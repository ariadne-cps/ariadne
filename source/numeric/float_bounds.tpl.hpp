/***************************************************************************
 *            float_bounds.tpl.hpp
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

template<class F> Nat Bounds<F>::output_places=8;

template<class F> Bounds<F>::Bounds(ExactDouble d, PR pr)
    : _l(d.get_d(),downward,pr),_u(d.get_d(),upward,pr) {
}

template<class F> Bounds<F>::Bounds(TwoExp t, PR pr)
    : _l(t,pr),_u(t,pr) {
}

template<class F> Bounds<F>::Bounds(Integer const& z, PR pr)
    : _l(z,downward,pr),_u(z,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Dyadic const& w, PR pr)
    : _l(w,downward,pr),_u(w,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Decimal const& d, PR pr)
    : Bounds(Rational(d),pr) {
}

template<class F> Bounds<F>::Bounds(Rational const& q, PR pr)
    : _l(q,downward,pr),_u(q,upward,pr) {
}

template<class F> Bounds<F>::Bounds(ExactDouble const& dl, ExactDouble const& du, PR pr)
    : _l(dl.get_d(),downward,pr),_u(du.get_d(),upward,pr) {
}

template<class F> Bounds<F>::Bounds(Dyadic const& wl, Dyadic const& wu, PR pr)
    : _l(wl,downward,pr),_u(wu,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Rational const& ql, Rational const& qu, PR pr)
    : _l(ql,downward,pr),_u(qu,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Bounds<F> const& x, PR pr)
    : _l(x._l,downward,pr), _u(x._u,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Real const& x, PR pr)
    : Bounds(x.get(pr)) {
}

template<class F> Bounds<F>::Bounds(LowerBound<F> const& lower, ValidatedUpperNumber const& upper)
    : Bounds<F>(lower,lower.create(upper)) { }

template<class F> Bounds<F>::Bounds(ValidatedLowerNumber const& lower, UpperBound<F> const& upper)
    : Bounds<F>(upper.create(lower),upper) { }

template<class F> Bounds<F>::Bounds(ValidatedLowerNumber const& lower, ValidatedUpperNumber const& upper, PR pr)
    : Bounds<F>(lower.get(LowerTag(),pr),upper.get(UpperTag(),pr)) { }

template<class F> Bounds<F>::Bounds(ValidatedNumber const& y, PR pr)
    : Bounds(y.get(BoundedTag(),pr)) {
}

template<class F> Bounds<F> Bounds<F>::create(ValidatedNumber const& y) const {
    return Bounds<F>(y,this->precision());
}

template<class F> Bounds<F>& Bounds<F>::operator=(ValidatedNumber const& y) {
    return *this = Bounds<F>(y,this->precision());
}

/*
template<class F> Bounds<F>::operator ValidatedNumber() const {
    return ValidatedNumber(new NumberWrapper<Bounds<F>>(*this));
}
*/

template<class F> Bounds<F> Operations<Bounds<F>>::_trunc(Bounds<F> const& x) {
    typename F::RoundingModeType rm=F::get_rounding_mode();
    const double& xl=x.lower_raw().get_d();
    const double& xu=x.upper_raw().get_d();
    // Use machine epsilon instead of minimum to move away from zero
    const float fm=std::numeric_limits<float>::epsilon();
    volatile float tu=xu;
    if(tu<xu) { F::set_rounding_upward(); tu+=fm; }
    volatile float tl=xl;
    if(tl>xl) { F::set_rounding_downward(); tl-=fm; }
    F::set_rounding_mode(rm);
    assert(tl<=xl); assert(tu>=xu);
    return Bounds<F>(double(tl),double(tu));
}

template<class F> Bounds<F> Operations<Bounds<F>>::_trunc(Bounds<F> const& x, Nat n) {
    Bounds<F> _e=Bounds<F>(std::pow(2.0,52-(Int)n));
    Bounds<F> y=x+_e; return y-_e;
}

template<class F> OutputStream& Operations<Bounds<F>>::_write(OutputStream& os, const Bounds<F>& x) {
    os << '{';
    write(os,x.lower().raw(),Bounds<F>::output_places,downward);
    os << ':';
    write(os,x.upper().raw(),Bounds<F>::output_places,upward);
    os << '}';
    return os;

}

template<class F> InputStream& Operations<Bounds<F>>::_read(InputStream& is, Bounds<F>& x) {
    char cl,cm,cr;
    F _l,_u;
    auto rnd=F::get_rounding_mode();
    is >> cl;
    F::set_rounding_downward();
    is >> _l;
    is >> cm;
    F::set_rounding_upward();
    is >> _u;
    is >> cr;
    F::set_rounding_mode(rnd);
    ARIADNE_ASSERT(not is.fail());
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    x._l=_l; x._u=_u;
    return is;
}



} // namespace Ariadne
