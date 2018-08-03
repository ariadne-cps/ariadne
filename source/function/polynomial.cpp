/***************************************************************************
 *            polynomial.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
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

template class AlgebraOperations<Polynomial<FloatDP>>;
template class AlgebraOperations<Polynomial<FloatDPApproximation>>;
template class AlgebraOperations<Polynomial<FloatDPBounds>>;
// template class AlgebraOperations<Polynomial<UpperIntervalType>>;

template<> Void Polynomial<FloatDPValue>::cleanup() { }

template Polynomial<FloatDPValue>::Polynomial(SizeType);
template Expansion<MultiIndex,FloatDPValue>& Polynomial<FloatDPValue>::expansion();
template OutputStream& Polynomial<FloatDPValue>::_write(OutputStream&) const;
template OutputStream& Polynomial<FloatDPValue>::_write(OutputStream&, List<String> const&) const;


template class Polynomial<FloatMPApproximation>;
template class Polynomial<FloatMPBounds>;
template class AlgebraOperations<Polynomial<FloatMPApproximation>>;
template class AlgebraOperations<Polynomial<FloatMPBounds>>;

} //namespace Ariadne


