/***************************************************************************
 *            expansion.cpp
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

#include "expansion.hpp"
#include "expansion.tpl.hpp"

#include "../numeric/float.hpp"
#include "../geometry/interval.hpp"


namespace Ariadne {
//    template class Expansion<DegreeType,FloatDP>;

    template class Expansion<MultiIndex,double>;
    template class Expansion<MultiIndex,Dyadic>;

    template class Expansion<MultiIndex,FloatDP>;
    template class Expansion<MultiIndex,FloatDPValue>;
    template class Expansion<MultiIndex,FloatDPBounds>;
    template class Expansion<MultiIndex,FloatDPApproximation>;

    template class SortedExpansion<MultiIndex,FloatDP,GradedIndexLess>;
    template class SortedExpansion<MultiIndex,FloatDPApproximation,GradedIndexLess>;
    template class SortedExpansion<MultiIndex,FloatDPBounds,GradedIndexLess>;

    template class SortedExpansion<MultiIndex,FloatDP,ReverseLexicographicIndexLess>;
    template class SortedExpansion<MultiIndex,FloatDPValue,ReverseLexicographicIndexLess>;
    template class SortedExpansion<MultiIndex,FloatDPBounds,ReverseLexicographicIndexLess>;
    template class SortedExpansion<MultiIndex,FloatDPApproximation,ReverseLexicographicIndexLess>;


    template class Expansion<MultiIndex,FloatDPUpperInterval>;
    template class SortedExpansion<MultiIndex,FloatDPUpperInterval,GradedIndexLess>;
    template class SortedExpansion<MultiIndex,FloatDPUpperInterval,ReverseLexicographicIndexLess>;


    template class Expansion<MultiIndex,FloatMP>;
    template class Expansion<MultiIndex,FloatMPValue>;
    template class Expansion<MultiIndex,FloatMPBounds>;
    template class Expansion<MultiIndex,FloatMPApproximation>;

    template class SortedExpansion<MultiIndex,FloatMP,GradedIndexLess>;
    template class SortedExpansion<MultiIndex,FloatMPApproximation,GradedIndexLess>;
    template class SortedExpansion<MultiIndex,FloatMPBounds,GradedIndexLess>;

    template class SortedExpansion<MultiIndex,FloatMP,ReverseLexicographicIndexLess>;
    template class SortedExpansion<MultiIndex,FloatMPValue,ReverseLexicographicIndexLess>;
    template class SortedExpansion<MultiIndex,FloatMPBounds,ReverseLexicographicIndexLess>;
    template class SortedExpansion<MultiIndex,FloatMPApproximation,ReverseLexicographicIndexLess>;

    template Expansion<MultiIndex,Interval<FloatDPUpperBound>>::Expansion(Expansion<MultiIndex,FloatDPBounds> const&);
}
