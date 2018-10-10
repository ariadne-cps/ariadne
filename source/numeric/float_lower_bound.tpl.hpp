/***************************************************************************
 *            float_lower_bound.tpl.hpp
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

template<class F> LowerBound<F>::LowerBound(ExactDouble d, PR pr)
    : _l(d.get_d(),pr) {
}

template<class F> LowerBound<F>::LowerBound(Integer const& z, PR pr)
    : _l(z,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(Dyadic const& w, PR pr)
    : _l(w,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(Decimal const& d, PR pr)
    : LowerBound(Rational(d),pr) {
}

template<class F> LowerBound<F>::LowerBound(Rational const& q, PR pr)
    : _l(q,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(Real const& r, PR pr)
    : LowerBound(r.get(pr)) {
}

template<class F> LowerBound<F>::LowerBound(LowerBound<F> const& x, PR pr)
    : _l(x._l,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(ValidatedLowerNumber const& y, PR pr)
    : LowerBound(y.get(LowerTag(),pr)) {
}

template<class F> LowerBound<F>& LowerBound<F>::operator=(ValidatedLowerNumber const& y) {
    return *this=LowerBound<F>(y,this->precision());
}

template<class F> LowerBound<F>::operator ValidatedLowerNumber() const {
    ARIADNE_NOT_IMPLEMENTED;
    //return ValidatedLowerNumber(new NumberWrapper<LowerBound<F>>(*this));
}

template<class F> LowerBound<F> LowerBound<F>::create(ValidatedLowerNumber const& y) const {
    return LowerBound<F>(y,this->precision());
}

template<class F> UpperBound<F> LowerBound<F>::create(ValidatedUpperNumber const& y) const {
    return UpperBound<F>(y,this->precision());
}



} // namespace Ariadne

