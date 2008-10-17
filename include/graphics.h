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
#include <string>

namespace Ariadne {

class Box;
class Polytope;

struct Colour {
  Colour();
  Colour(unsigned char rd, unsigned char gr, unsigned char bl, bool tr=true);
  Colour(const char* nm, unsigned char rd, unsigned char gr, unsigned char bl, bool tr=true);
  std::string name;
  unsigned char red, green, blue;
  bool transparant;
};

std::ostream& operator<<(std::ostream& os, const Colour& c);

extern const Colour transparant;

extern const Colour white;
extern const Colour black;
extern const Colour red;
extern const Colour green;
extern const Colour blue;
extern const Colour yellow;
extern const Colour cyan;
extern const Colour magenta;

struct LineStyle { explicit LineStyle(bool ls) : _style(ls) { } operator bool() const { return this->_style; } private: bool _style; };
struct LineWidth { explicit LineWidth(double lw) : _width(lw) { } operator double() const { return this->_width; } private: double _width; };
struct LineColour : Colour { LineColour(const Colour& lc) : Colour(lc) { } };
struct FillColour : Colour { FillColour(const Colour& fc) : Colour(fc) { } };

inline LineStyle line_style(bool s) { return LineStyle(s); }
inline LineWidth line_width(double w) { return LineWidth(w); }
inline LineColour line_colour(const Colour& c) { return LineColour(c); }
inline FillColour fill_colour(const Colour& c) { return FillColour(c); }
    

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
