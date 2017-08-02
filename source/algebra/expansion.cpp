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

    template class Expansion<FloatDP>;
    template class Expansion<FloatDPValue>;
    template class Expansion<FloatDPBounds>;
    template class Expansion<FloatDPApproximation>;

    template class SortedExpansion<FloatDP,GradedIndexLess>;
    template class SortedExpansion<FloatDPApproximation,GradedIndexLess>;
    template class SortedExpansion<FloatDPBounds,GradedIndexLess>;

    template class SortedExpansion<FloatDP,ReverseLexicographicIndexLess>;
    template class SortedExpansion<FloatDPValue,ReverseLexicographicIndexLess>;
    template class SortedExpansion<FloatDPBounds,ReverseLexicographicIndexLess>;
    template class SortedExpansion<FloatDPApproximation,ReverseLexicographicIndexLess>;


    template class Expansion<FloatDPUpperInterval>;
    template class SortedExpansion<FloatDPUpperInterval,GradedIndexLess>;
    template class SortedExpansion<FloatDPUpperInterval,ReverseLexicographicIndexLess>;


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
