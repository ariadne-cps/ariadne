/***************************************************************************
 *            expansion.cpp
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
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "expansion.hpp"
#include "expansion.tpl.hpp"

#include "numeric/float.hpp"
#include "geometry/interval.hpp"


namespace Ariadne {
    template class Expansion<double>;
    template class Expansion<Dyadic>;

    template class Expansion<Float64>;
    template class Expansion<Float64Value>;
    template class Expansion<Float64Bounds>;
    template class Expansion<Float64Approximation>;

    template class SortedExpansion<Float64,GradedIndexLess>;
    template class SortedExpansion<Float64Approximation,GradedIndexLess>;
    template class SortedExpansion<Float64Bounds,GradedIndexLess>;

    template class SortedExpansion<Float64,ReverseLexicographicIndexLess>;
    template class SortedExpansion<Float64Value,ReverseLexicographicIndexLess>;
    template class SortedExpansion<Float64Bounds,ReverseLexicographicIndexLess>;
    template class SortedExpansion<Float64Approximation,ReverseLexicographicIndexLess>;


    template class Expansion<Float64UpperInterval>;
    template class SortedExpansion<Float64UpperInterval,GradedIndexLess>;
    template class SortedExpansion<Float64UpperInterval,ReverseLexicographicIndexLess>;


    template class Expansion<FloatMP>;
    template class Expansion<FloatMPValue>;
    template class Expansion<FloatMPBounds>;
    template class Expansion<FloatMPApproximation>;

    template class SortedExpansion<FloatMP,GradedIndexLess>;
    template class SortedExpansion<FloatMPApproximation,GradedIndexLess>;
    template class SortedExpansion<FloatMPBounds,GradedIndexLess>;

    template class SortedExpansion<FloatMP,ReverseLexicographicIndexLess>;
    template class SortedExpansion<FloatMPValue,ReverseLexicographicIndexLess>;
    template class SortedExpansion<FloatMPBounds,ReverseLexicographicIndexLess>;
    template class SortedExpansion<FloatMPApproximation,ReverseLexicographicIndexLess>;

}
