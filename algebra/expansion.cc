/***************************************************************************
 *            expansion.cc
 *
 *  Copyright 2008-15  Pieter Collins
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

#include "expansion.h"
#include "expansion.tcc"

#include "numeric/float.h"
#include "geometry/interval.h"

namespace Ariadne {
    template class Expansion<Float>;
    template class Expansion<ExactFloat>;
    template class Expansion<ValidatedFloat>;
    template class Expansion<ApproximateFloat>;
    template class Expansion<UpperInterval>;
    template class Expansion<ExactInterval>;

    template class SortedExpansion<Float,GradedKeyLess>;
    template class SortedExpansion<ApproximateFloat,GradedKeyLess>;
    template class SortedExpansion<ValidatedFloat,GradedKeyLess>;
    template class SortedExpansion<UpperInterval,GradedKeyLess>;

    template class SortedExpansion<Float,ReverseLexicographicKeyLess>;
    template class SortedExpansion<ExactFloat,ReverseLexicographicKeyLess>;
    template class SortedExpansion<ValidatedFloat,ReverseLexicographicKeyLess>;
    template class SortedExpansion<ApproximateFloat,ReverseLexicographicKeyLess>;
    template class SortedExpansion<UpperInterval,ReverseLexicographicKeyLess>;
    template class SortedExpansion<ExactInterval,ReverseLexicographicKeyLess>;
}

#include "algebra/differential.h"

namespace Ariadne {
}
