/***************************************************************************
 *            polynomial.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "../numeric/numeric.hpp"
#include "../config.hpp"

#include "../geometry/interval.hpp"
#include "../function/polynomial.hpp"
#include "../function/polynomial.tpl.hpp"

namespace Ariadne {

template class Polynomial<FloatDP>;
template class Polynomial<FloatDPApproximation>;
template class Polynomial<FloatDPBounds>;
template class Polynomial<UpperIntervalType>;

template struct AlgebraOperations<Polynomial<FloatDP>>;
template struct AlgebraOperations<Polynomial<FloatDPApproximation>>;
template struct AlgebraOperations<Polynomial<FloatDPBounds>>;
// template struct AlgebraOperations<Polynomial<UpperIntervalType>>;

template<> Void Polynomial<FloatDPValue>::cleanup() { }

template Polynomial<FloatDPValue>::Polynomial(SizeType);
template Expansion<MultiIndex,FloatDPValue>& Polynomial<FloatDPValue>::expansion();
template OutputStream& Polynomial<FloatDPValue>::_write(OutputStream&) const;
template OutputStream& Polynomial<FloatDPValue>::_write(OutputStream&, List<String> const&) const;


template class Polynomial<FloatMPApproximation>;
template class Polynomial<FloatMPBounds>;
template struct AlgebraOperations<Polynomial<FloatMPApproximation>>;
template struct AlgebraOperations<Polynomial<FloatMPBounds>>;

} //namespace Ariadne


