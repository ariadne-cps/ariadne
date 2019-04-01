/***************************************************************************
 *            chebyshev_model.cpp
 *
 *  Copyright 2008--18  Pieter Collins
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
#include "algebra/vector.hpp"

#include "function/chebyshev_model.hpp"
#include "function/chebyshev_model.tpl.hpp"

#include "function/scaled_function_patch.hpp"
#include "function/scaled_function_patch.tpl.hpp"

namespace Ariadne {

template class MultivariateChebyshevModel<FloatDPApproximation>;
template class MultivariateChebyshevModel<FloatMPApproximation>;

template class ScaledFunctionPatch<MultivariateChebyshevModel<FloatDPApproximation>>;
template class VectorScaledFunctionPatch<MultivariateChebyshevModel<FloatDPApproximation>>;
template class ScaledFunctionPatch<MultivariateChebyshevModel<FloatMPApproximation>>;
template class VectorScaledFunctionPatch<MultivariateChebyshevModel<FloatMPApproximation>>;

} // namespace Ariadne
