/***************************************************************************
 *            output/colour.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file output/colour.hpp
 *  \brief Colours for graphics objects.
 */

#ifndef ARIADNE_COLOUR_HPP
#define ARIADNE_COLOUR_HPP

#include <iosfwd>
#include <string>

typedef unsigned int Nat;

namespace Ariadne {


struct Colour {
    Colour();
    Colour(double rd, double gr, double bl, double op=1.0);
    Colour(double rd, double gr, double bl, Bool tr);
    Colour(const char* nm, double rd, double gr, double bl, double op=1.0);
    Colour(const char* nm, double rd, double gr, double bl, Bool tr);
    StringType name;
    double red, green, blue;
    double opacity;
};

OutputStream& operator<<(OutputStream& os, const Colour& c);

extern const Colour transparant;

extern const Colour white;
extern const Colour black;
extern const Colour red;
extern const Colour green;
extern const Colour blue;
extern const Colour yellow;
extern const Colour cyan;
extern const Colour magenta;
extern const Colour ariadneorange;


} // namespace Ariadne

#endif // ARIADNE_COLOUR_HPP
