/***************************************************************************
 *            numeric_examples.cpp
 *
 *  Copyright  2020  Pieter Collins
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

#include <ariadne.hpp>
#include "ariadne/ariadne.hpp"

using namespace Ariadne;


int main() {

    //! [Real usage]
    Real r=6*asin(sqrt(rec(3));
    FloatMPBounds x=r.compute(Accuracy(128_bits)).get(precision(96_bits));
    //! [Real usage]

}
