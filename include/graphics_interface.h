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

class DrawableInterface;
class FigureInterface;
class CanvasInterface;
class PlanarProjectionMap;


struct Vector2d { double x,y; Vector2d(double xx, double yy) : x(xx), y(yy) { } };
inline Vector2d operator-(const Vector2d& v) { return Vector2d(-v.x,-v.y); }
inline Vector2d operator+(const Vector2d& v1, const Vector2d& v2) { return Vector2d(v1.x+v2.x,v1.y+v2.y); }
inline Vector2d operator-(const Vector2d& v1, const Vector2d& v2) { return Vector2d(v1.x-v2.x,v1.y-v2.y); }
inline Vector2d operator*(const double& s1, const Vector2d& v2) { return Vector2d(s1*v2.x,s1*v2.y); }

struct Point2d { double x,y; Point2d(double xx, double yy) : x(xx), y(yy) { } };
inline Point2d& operator+=(Point2d& pt, const Vector2d& v) { pt.x+=v.x; pt.y+=v.y; return pt; }
inline Point2d& operator-=(Point2d& pt, const Vector2d& v) { pt.x-=v.x; pt.y-=v.y; return pt; }

//! \brief Base interface for plotting and drawing classes.
class FigureInterface {
  public:
    virtual ~FigureInterface() { };
    virtual void set_projection(uint as, uint ix, uint iy) { };
    virtual void set_x_axis_label(const string&) { };
    virtual void set_y_axis_label(const string&) { };
    virtual void set_line_style(bool) { };
    virtual void set_line_width(double) { };
    virtual void set_line_colour(Colour) { };
    virtual void set_fill_colour(Colour) { };
    virtual string get_x_axis_label() const { return ""; };
    virtual string get_y_axis_label() const { return ""; };
    virtual bool get_line_style() const { return true; };
    virtual double get_line_width() const { return 1.0; };
    virtual Colour get_line_colour() const { return black; };
    virtual bool get_fill_style() const { return true; };
    virtual Colour get_fill_colour() const { return white; };
    virtual void draw(const DrawableInterface&) = 0;
};

inline void draw(FigureInterface& fig, const DrawableInterface& shape) {
    fig.draw(shape);
}

class CanvasInterface {
  public:
    virtual ~CanvasInterface() { };
    virtual uint x_coordinate() const = 0;
    virtual uint y_coordinate() const = 0;
    virtual void move_to(double x, double y) = 0;
    virtual void line_to(double x, double y) = 0;
    virtual void circle(double x, double y, double r) = 0;
    virtual void dot(double x, double y) = 0;
    virtual void stroke() = 0;
    virtual void fill() = 0;
    virtual void clip() = 0;
    virtual double get_line_width() const = 0;
    virtual void set_line_width(double lw) = 0;
    virtual void set_line_colour(double r, double g, double b) = 0;
    virtual void set_fill_colour(double r, double g, double b) = 0;
    virtual void set_bounding_box(double x0, double x1, double y0, double y1) = 0;
};


//! \brief Base interface for drawable objects
class DrawableInterface {
  public:
    virtual ~DrawableInterface() { }
    virtual DrawableInterface* clone() const = 0;
    virtual void draw(CanvasInterface& c) const = 0;
    virtual uint dimension() const = 0;
    virtual std::ostream& write(std::ostream& os) const { return os << "Drawable"; }
};


} // namespace Ariadne


#endif // ARIADNE_GRAPHICS_INTERFACE_H
