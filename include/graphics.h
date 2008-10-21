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

#include "box.h"
#include "zonotope.h"
#include "approximate_taylor_model.h"
#include "geometry2d.h"

namespace Ariadne {

class Box;
class Polytope;
class PlanarProjectionMap;

struct Colour {
  Colour();
  Colour(double rd, double gr, double bl, bool tr=true);
  Colour(const char* nm, double rd, double gr, double bl, bool tr=true);
  std::string name;
  double red, green, blue;
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
struct FillStyle { explicit FillStyle(bool fs) : _style(fs) { } operator bool() const { return this->_style; } private: bool _style; };
struct FillColour : Colour { FillColour(const Colour& fc) : Colour(fc) { } };

inline LineStyle line_style(bool s) { return LineStyle(s); }
inline LineWidth line_width(double w) { return LineWidth(w); }
inline LineColour line_colour(const Colour& c) { return LineColour(c); }
inline FillStyle fill_style(bool s) { return FillStyle(s); }
inline FillColour fill_colour(const Colour& c) { return FillColour(c); }
    

class Graphic {
  public:
    ~Graphic();
    Graphic();
    void set_projection_map(const PlanarProjectionMap&);
    void set_bounding_box(const Box&);
    void set_line_style(bool);
    void set_line_width(double);
    void set_line_colour(Colour);
    void set_fill_style(bool);
    void set_fill_colour(Colour);
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


Graphic& operator<<(Graphic& g, const LineStyle& ls) { g.set_line_style(ls); return g; }
Graphic& operator<<(Graphic& g, const LineWidth& lw) { g.set_line_width(lw); return g; }
Graphic& operator<<(Graphic& g, const LineColour& lc) { g.set_line_colour(lc); return g; }
Graphic& operator<<(Graphic& g, const FillStyle& fs) { g.set_fill_style(fs); return g; }
Graphic& operator<<(Graphic& g, const FillColour& fc) { g.set_fill_colour(fc); return g; }

inline void plot(Graphic& g, const Box& bx) { g.plot(bx); }
inline void plot(Graphic& g, const Zonotope& z) { g.plot(polytope(z)); }
inline void plot(Graphic& g, const Polytope& p) { g.plot(polytope(p)); }
inline void plot(Graphic& g, const ApproximateTaylorModel& ts) { g.plot(polytope(ts)); }
template<class BS> inline void plot(Graphic& g, const std::pair<int,BS>& hs) { plot(g,hs.second); }

template<class SET> Graphic& operator<<(Graphic& g, const SET& set) { plot(g,set); return g; }




} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
