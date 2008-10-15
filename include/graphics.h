/***************************************************************************
 *            graphics.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file graphics.h
 *  \brief Graphics class for plotting and outputting figures.
 */

#ifndef ARIADNE_GRAPHICS_H
#define ARIADNE_GRAPHICS_H

#include <iosfwd>

namespace Ariadne {

class Box;
class Polytope;

enum Colour { transparant,black,white,red,green,blue,yellow,cyan,magenta };
Colour fill_colour(Colour c) { return c; }
double line_width(double l) { return l; }
bool line_style(bool l) { return l; }

class Graphic {
  public:
    ~Graphic();
    Graphic();
    Graphic(std::ofstream& ofs);
    void set_bounding_box(const Box&);
    void plot(const Box&);
    void plot(const Polytope&);
    void clear();
    void display();
    void write(const char* filename);
  public:
    class Impl;
  private:
    Impl* _impl;
};


template<class X> Graphic& operator<<(Graphic& g, const X& x) { return g; }
template<> Graphic& operator<<(Graphic& g, const Box& bx);

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
