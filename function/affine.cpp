/***************************************************************************
 *            affine.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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

#include "numeric/numeric.hpp"
#include "config.h"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/affine.hpp"

namespace Ariadne {

template class Affine<Float64Approximation>;
template class Affine<Float64Bounds>;
template struct ProvideAlgebraOperations<Affine<Float64Approximation>,Float64Approximation>;
template struct ProvideAlgebraOperations<Affine<Float64Bounds>,Float64Bounds>;

template class Affine<FloatMPApproximation>;
template class Affine<FloatMPBounds>;
template struct ProvideAlgebraOperations<Affine<FloatMPApproximation>,FloatMPApproximation>;
template struct ProvideAlgebraOperations<Affine<FloatMPBounds>,FloatMPBounds>;

} //namespace Ariadne


