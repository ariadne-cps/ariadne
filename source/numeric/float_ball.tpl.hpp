/***************************************************************************
 *            float_ball.tpl.hpp
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

template<class F, class FE> Ball<F,FE>::Ball(ExactDouble d, PR pr)
    : _v(d.get_d(),to_nearest,pr), _e(0,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(TwoExp t, PR pr)
    : _v(t,pr), _e(0u,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Integer const& z, PR pr) : Ball<F,FE>(Rational(z),pr) {
}

template<class F, class FE> Ball<F,FE>::Ball(Dyadic const& w, PR pr)
    : _v(F(w,to_nearest,pr)), _e(abs(Dyadic(_v)-w),upward,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Decimal const& d, PR pr)
    : Ball(Rational(d),pr) {
}

template<class F, class FE> Ball<F,FE>::Ball(Rational const& q, PR pr)
    : _v(F(q,to_nearest,pr)), _e(abs(Rational(_v)-q),upward,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Real const& r, PR pr)
    : Ball(r.get(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Ball<F,FE> const& x, PR pr)
    : _v(x._v,to_nearest,pr), _e(x._e,upward,_error_precision<PRE>(pr))
{
    F d = (this->_v>=x._v) ? sub(up,this->_v,x._v) : sub(up,x._v,this->_v);
    _e=add(up,_e,_make_error<FE>(d));
}

template<class F, class FE> Ball<F,FE>::Ball(ValidatedNumber const& y, PR pr)
    : Ball(y.get(MetricTag(),pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Bounds<F> const& x, PRE pre)
    : _v(x.value_raw()) , _e(_make_error<FE>(x.error_raw(),pre)) {
}



template<class F, class FE> Ball<F,FE>::Ball(Real const& r, PR pr, PRE pre)
    : Ball(r.get(pr),pre) {
}

template<class F, class FE> Ball<F,FE>::Ball(Rational const& q, PR pr, PRE pre)
    : _v(F(q,to_nearest,pr)), _e(abs(Rational(_v)-q),upward,pre) {
}

template<class F, class FE> Ball<F,FE>::Ball(ValidatedNumber const& y, PR pr, PRE pre)
    : Ball(y.get(MetricTag(),pr,pre)) {
}

template<class F, class FE> Ball<F,FE> Ball<F,FE>::create(ValidatedNumber const& y) const {
    return Ball<F,FE>(y,this->precision());
}

template<class F, class FE> Ball<F,FE>& Ball<F,FE>::operator=(ValidatedNumber const& y)
{
    return *this = Ball(y,this->precision(),this->error_precision());
}

/*
template<class F, class FE> Ball<F,FE>::operator ValidatedNumber() const {
    return ValidatedNumber(new NumberWrapper<Ball<F,FE>>(*this));
}
*/



template<class F, class FE> OutputStream& Operations<Ball<F,FE>>::_write(OutputStream& os, Ball<F,FE> const& x) {
    return os << x.value() << "\u00b1" << x.error();
}

template<class F, class FE> InputStream& Operations<Ball<F,FE>>::_read(InputStream& is, Ball<F,FE>& x) {
    static const char pmstr[] = "\u00b1";
    char cpm[3];
    F _v; FE _e;
    auto rnd=F::get_rounding_mode();
    F::set_rounding_to_nearest();
    is >> _v;
    is >> cpm[0] >> cpm[1];
    F::set_rounding_mode(rnd);
    auto rnde=FE::get_rounding_mode();
    FE::set_rounding_upward();
    is >> _e;
    FE::set_rounding_mode(rnde);
    ARIADNE_ASSERT(not is.fail());
    ARIADNE_ASSERT(std::strcmp(cpm,pmstr));
    x._v=_v; x._e=_e;
    return is;
}



} // namespace Ariadne
