/***************************************************************************
 *            regular_expression.cpp
 *
 *  Copyright  2018  Pieter Collins
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

#include "numeric/integer.hpp"
#include "regular_expression.tpl.hpp"

namespace Ariadne {

template RegularExpression<Char> simplify(RegularExpression<Char>);
template RegularExpression<Pair<Char,Char>> simplify(RegularExpression<Pair<Char,Char>>);

template RegularExpression<Int16> simplify(RegularExpression<Int16>);
template RegularExpression<Pair<Int16,Int16>> simplify(RegularExpression<Pair<Int16,Int16>>);


} // namespace Ariadne

