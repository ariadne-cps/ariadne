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

const FloatDPValue infty = FloatDPValue(FloatDP::inf(dp));

FloatValue<DoublePrecision> operator"" _exact(long double lx) {
    double x=lx;
    assert(x==lx);
    return FloatValue<DoublePrecision>(x);
}

template class Value<FloatDP>;
template class Operations<Value<FloatDP>>;
template class Value<FloatMP>;
template class Operations<Value<FloatMP>>;

template<> String class_name<Value<FloatDP>>() { return "FloatDPValue"; }
template<> String class_name<Value<FloatMP>>() { return "FloatMPValue"; }


} // namespace Ariadne
