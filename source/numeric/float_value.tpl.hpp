/***************************************************************************
 *            float_value.tpl.hpp
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


#ifndef ARIADNE_FLOAT_VALUE_TPL_HPP
#define ARIADNE_FLOAT_VALUE_TPL_HPP

#include "float_value.hpp"

#include "float_ball.hpp"

#include "builtin.hpp"

namespace Ariadne {

template<class F> Nat Value<F>::output_places = 17;

template<class F> Value<F>::Value(ExactDouble const& d, PR pr)
    : _v(d,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==d,"Exact double"<<d<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(TwoExp const& t, PR pr)
    : _v(t,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==Dyadic(t),"Power-of-two "<<t<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Integer const& z, PR pr)
    : _v(z,to_nearest,pr)
{
    Rational q(_v);
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==z,"Integer "<<z<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Dyadic const& w, PR pr)
    : _v(w,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==w || is_nan(w),"Dyadic number "<<w<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Value<F> const& x, PR pr)
    : _v(x._v,to_nearest,pr)
{
    ARIADNE_ASSERT_MSG(*this==x,"Exact FloatValue "<<x<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::operator Dyadic() const {
    return this->_v.operator Dyadic();
}

template<class F> Value<F>::operator Rational() const {
    return Rational(this->operator Dyadic());
}

template<class F> Value<F>& Value<F>::operator=(TwoExp const& t) {
    _v=F(t,this->precision());
    return *this;
}

template<class F> Value<F>& Value<F>::operator=(Integer const& z) {
    _v=F(z,this->precision());
    ARIADNE_ASSERT_MSG(Dyadic(_v)==z,"Integer "<<z<<" cannot be assigned exactly to a floating-point number with precision "<<this->precision()<<"; nearest is "<<(*this));
    return *this;
}

template<class F> Value<F>& Value<F>::operator=(Dyadic const& w) {
    _v=F(w,this->precision());
    ARIADNE_ASSERT_MSG(Dyadic(_v)==w,"Dyadic number "<<w<<" cannot be assigned exactly to a floating-point number with precision "<<this->precision()<<"; nearest is "<<(*this));
    return *this;
}

template<class F> auto Value<F>::create(ValidatedNumber const& y) const -> Ball<F> {
    return Ball<F>(y,this->precision());
}

/*
template<class F> Value<F>::operator ExactNumber() const {
    return ExactNumber(new NumberWrapper<Value<F>>(*this));
}
*/

} // namespace Ariadne

#endif
