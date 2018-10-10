/***************************************************************************
 *            float_approximation.tpl.hpp
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

template<class F> Nat Approximation<F>::output_places = 4;

template<class F> Approximation<F>::Approximation(double d, PR pr)
    : _a(d,to_nearest,pr)
{
}

template<class F> Approximation<F>::Approximation(ExactDouble d, PR pr)
    : _a(d.get_d(),to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(TwoExp t, PR pr)
    : _a(t,pr) {
}

template<class F> Approximation<F>::Approximation(Integer const& z, PR pr)
    : _a(z,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Dyadic const& w, PR pr)
    : _a(w,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Decimal const& d, PR pr)
    : Approximation<F>(Rational(d),pr) {
}

template<class F> Approximation<F>::Approximation(Rational const& q, PR pr)
    : _a(q,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Approximation<F> const& x, PR pr)
    : _a(x._a,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Real const& r, PR pr)
    : Approximation<F>(r.get(pr)) {
}

template<class F> Approximation<F>::Approximation(ApproximateNumber const& y, PR pr)
    : Approximation<F>(y.get(ApproximateTag(),pr)) {
}

template<class F> Approximation<F>& Approximation<F>::operator=(ApproximateNumber const& y) {
    return *this=Approximation<F>(y,this->precision());
}

template<class F> Approximation<F> Approximation<F>::create(ApproximateNumber const& y) const {
    return Approximation<F>(y,this->precision());
}
/*
template<class F> Approximation<F>::operator ApproximateNumber() const {
    return ApproximateNumber(new NumberWrapper<Approximation<F>>(*this));
}
*/



} // namespace Ariadne
