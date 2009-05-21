/***************************************************************************
 *            graphics_interface.h
 *
 *  Copyright 2009  Davide Bresolin
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
 
/*! \file graphics_interface.h
 *  \brief Base graphics interface from which all plotting and drawing classes are inherited.
 */

#ifndef ARIADNE_GRAPHICS_INTERFACE_H
#define ARIADNE_GRAPHICS_INTERFACE_H

#include "macros.h"
#include "colour.h"

typedef unsigned int uint;

namespace Ariadne {

using namespace std;

class Point;
class Box;
class Polytope;
class InterpolatedCurve;
class Zonotope;

class Colour;

//! \brief Base interface for plotting and drawing classes.
class GraphicsInterface {
  public:
    virtual ~GraphicsInterface() { };
    virtual void set_line_style(bool) { };
    virtual void set_line_width(double) { };
    virtual void set_line_colour(Colour) { };
    virtual void set_fill_colour(Colour) { };
    virtual bool get_line_style() const { return true; };
    virtual double get_line_width() const { return 1.0; };
    virtual Colour get_line_colour() const { return black; };
    virtual bool get_fill_style() const { return true; };
    virtual Colour get_fill_colour() const { return white; };
    virtual void draw(const std::vector<Point>&) = 0; // Draw a shape bounded by a list of points
    virtual void draw(const Point&) = 0;
    virtual void draw(const Box&) = 0;
    virtual void draw(const Polytope&) = 0;
    virtual void draw(const InterpolatedCurve&) = 0;
};

inline void draw(GraphicsInterface& g, const Point& pt) { g.draw(pt); }
inline void draw(GraphicsInterface& g, const Box& bx) { g.draw(bx); }
inline void draw(GraphicsInterface& g, const Polytope& p) { g.draw(p); }
inline void draw(GraphicsInterface& g, const InterpolatedCurve& c) { g.draw(c); }

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_INTERFACE_H
