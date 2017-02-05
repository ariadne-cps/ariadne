/***************************************************************************
 *            polynomial.cc
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

#include "numeric/numeric.h"
#include "config.h"

#include "geometry/interval.h"
#include "function/polynomial.h"
#include "function/polynomial.tpl.h"

namespace Ariadne {

template class Polynomial<Float64>;
template class Polynomial<Float64Approximation>;
template class Polynomial<Float64Bounds>;
template class Polynomial<UpperIntervalType>;

template class AlgebraOperations<Polynomial<Float64>>;
template class AlgebraOperations<Polynomial<Float64Approximation>>;
template class AlgebraOperations<Polynomial<Float64Bounds>>;
// template class AlgebraOperations<Polynomial<UpperIntervalType>>;

template<> Void Polynomial<Float64Value>::cleanup() { }

template Polynomial<Float64Value>::Polynomial(SizeType);
template Expansion<Float64Value>& Polynomial<Float64Value>::expansion();
template OutputStream& Polynomial<Float64Value>::_write(OutputStream&) const;
template OutputStream& Polynomial<Float64Value>::_write(OutputStream&, List<String> const&) const;


template class Polynomial<FloatMPApproximation>;
template class Polynomial<FloatMPBounds>;
template class AlgebraOperations<Polynomial<FloatMPApproximation>>;
template class AlgebraOperations<Polynomial<FloatMPBounds>>;

} //namespace Ariadne


