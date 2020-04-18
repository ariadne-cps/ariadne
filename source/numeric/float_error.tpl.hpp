/***************************************************************************
 *            float_error.tpl.hpp
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


#ifndef ARIADNE_FLOAT_ERROR_TPL_HPP
#define ARIADNE_FLOAT_ERROR_TPL_HPP

#include "float_error.hpp"

#include "float_value.hpp"
#include "float_bounds.hpp"

namespace Ariadne {

template<class F> Nat Error<F>::output_places = 3;

template<class F> Error<F>::Error(PositiveBounds<F> const& x)
    : _e(x._u) { }

template<class F> Error<F>::Error(Positive<Value<F>> const& x)
    : _e(x._v) { }

template<class F> struct Operations<Error<F>> {
    static OutputStream& _write(OutputStream& os, Error<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Error<F>::output_places},upward);
    }
    static InputStream& _read(InputStream& is, Error<F>& x) {
        UpperBound<F> xu; is >> xu; x=Error<F>(xu); return is;
    }
};


} // namespace Ariadne

#endif
