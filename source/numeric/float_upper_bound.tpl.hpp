/***************************************************************************
 *            float_upper_bound.tpl.hpp
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

template<class F> UpperBound<F>::UpperBound(ExactDouble d, PR pr)
    : _u(d.get_d(),upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(TwoExp t, PR pr)
    : _u(t,pr) {
}

template<class F> UpperBound<F>::UpperBound(Integer const& z, PR pr)
    : _u(z,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(Dyadic const& w, PR pr)
    : _u(w,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(Decimal const& d, PR pr)
    : UpperBound(Rational(d),pr) {
}

template<class F> UpperBound<F>::UpperBound(Rational const& q, PR pr)
    : _u(q,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(Real const& r, PR pr)
    : UpperBound(r.get(pr)) {
}

template<class F> UpperBound<F>::UpperBound(UpperBound<F> const& x, PR pr)
    : _u(x._u,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(ValidatedUpperNumber const& y, PR pr)
    : UpperBound(y.get(UpperTag(),pr)) {
}

template<class F> UpperBound<F>& UpperBound<F>::operator=(ValidatedUpperNumber const& y) {
    return *this = UpperBound<F>(y,this->precision());
}

template<class F> UpperBound<F>::operator ValidatedUpperNumber() const {
    ARIADNE_NOT_IMPLEMENTED;
    // return ValidatedUpperNumber(new NumberWrapper<UpperBound<F>>(*this));
}

template<class F> LowerBound<F> UpperBound<F>::create(ValidatedLowerNumber const& y) const {
    return LowerBound<F>(y,this->precision());
}

template<class F> UpperBound<F> UpperBound<F>::create(ValidatedUpperNumber const& y) const {
    return UpperBound<F>(y,this->precision());
}



} // namespace Ariadne
