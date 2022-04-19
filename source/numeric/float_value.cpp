/***************************************************************************
 *            numeric/float_value.cpp
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

#include "float_value.hpp"
#include "float_value.tpl.hpp"

#include "floatdp.hpp"
#include "floatmp.hpp"

namespace Ariadne {

const Dyadic infty = Dyadic::inf();

template<> Bool same(FloatDP const& x1, FloatDP const& x2) { return x1==x2; }
template<> Bool same(FloatMP const& x1, FloatMP const& x2) { return x1==x2; }

inline FloatDP const& cast_exact(UpperBound<FloatDP> const& x) { return reinterpret_cast<FloatDP const&>(x); }
inline FloatMP const& cast_exact(UpperBound<FloatMP> const& x) { return reinterpret_cast<FloatMP const&>(x); }
inline FloatDP const& cast_exact(Approximation<FloatDP> const& x) { return reinterpret_cast<FloatDP const&>(x); }
inline FloatMP const& cast_exact(Approximation<FloatMP> const& x) { return reinterpret_cast<FloatMP const&>(x); }

inline UpperBound<FloatDP> const& cast_unsigned(Positive<UpperBound<FloatDP>> const& x) { return reinterpret_cast<UpperBound<FloatDP> const&>(x); }
inline UpperBound<FloatMP> const& cast_unsigned(Positive<UpperBound<FloatMP>> const& x) { return reinterpret_cast<UpperBound<FloatMP> const&>(x); }
inline Approximation<FloatDP> const& cast_unsigned(Positive<Approximation<FloatDP>> const& x) { return reinterpret_cast<Approximation<FloatDP> const&>(x); }
inline Approximation<FloatMP> const& cast_unsigned(Positive<Approximation<FloatMP>> const& x) { return reinterpret_cast<Approximation<FloatMP> const&>(x); }

template<> Positive<FloatDP> cast_exact(Positive<UpperBound<FloatDP>> const& x) {
    return cast_positive(cast_exact(cast_unsigned(x))); }
template<> Positive<FloatMP> cast_exact(Positive<UpperBound<FloatMP>> const& x) {
    return cast_positive(cast_exact(cast_unsigned(x))); }

template<> Positive<FloatDP> cast_exact(Positive<Approximation<FloatDP>> const& x) {
    return cast_positive(cast_exact(cast_unsigned(x))); }
template<> Positive<FloatMP> cast_exact(Positive<Approximation<FloatMP>> const& x) {
    return cast_positive(cast_exact(cast_unsigned(x))); }

Bounds<FloatMP> add(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return add(Bounds<FloatMP>(x1),Bounds<FloatMP>(x2)); }
Bounds<FloatMP> sub(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return sub(Bounds<FloatMP>(x1),Bounds<FloatMP>(x2)); }
Bounds<FloatMP> mul(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return mul(Bounds<FloatMP>(x1),Bounds<FloatMP>(x2)); }
Bounds<FloatMP> div(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return div(Bounds<FloatMP>(x1),Bounds<FloatMP>(x2)); }

} // namespace Ariadne
