/***************************************************************************
 *      graphics_interface.h
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

template<class R, class A> inline R numeric_cast(const A&);

struct Vector2d {
    double x,y; Vector2d(double xx, double yy) : x(xx), y(yy) { }
    template<class X, class Y> Vector2d(const X& xx, const Y& yy) : x(numeric_cast<double>(xx)), y(numeric_cast<double>(yy)) { }
};
inline Vector2d operator-(const Vector2d& v) { return Vector2d(-v.x,-v.y); }
inline Vector2d operator+(const Vector2d& v1, const Vector2d& v2) { return Vector2d(v1.x+v2.x,v1.y+v2.y); }
inline Vector2d operator-(const Vector2d& v1, const Vector2d& v2) { return Vector2d(v1.x-v2.x,v1.y-v2.y); }
inline Vector2d operator*(const double& s1, const Vector2d& v2) { return Vector2d(s1*v2.x,s1*v2.y); }
inline std::ostream& operator<<(std::ostream& os, const Vector2d& v) { return os << "["<<v.x<<","<<v.y<<"]"; }

struct Point2d {
    double x,y; Point2d(double xx, double yy) : x(xx), y(yy) { }
    template<class X, class Y> Point2d(const X& xx, const Y& yy) : x(numeric_cast<double>(xx)), y(numeric_cast<double>(yy)) { }
};
inline bool operator==(Point2d& pt1, const Point2d& pt2) { return pt1.x==pt2.x && pt1.y==pt2.y; }
inline Point2d& operator+=(Point2d& pt, const Vector2d& v) { pt.x+=v.x; pt.y+=v.y; return pt; }
inline Point2d& operator-=(Point2d& pt, const Vector2d& v) { pt.x-=v.x; pt.y-=v.y; return pt; }
inline std::ostream& operator<<(std::ostream& os, const Point2d& pt) { return os << "("<<pt.x<<","<<pt.y<<")"; }

struct Box2d { double xl,xu,yl,yu; Box2d() { } Box2d(double xxl, double xxu, double yyl, double yyu) : xl(xxl), xu(xxu), yl(yl), yu(yyu) { } };
inline std::ostream& operator<<(std::ostream& os, const Box2d& bx) { return os << "["<<bx.xl<<","<<bx.xu<<"]x["<<bx.yl<<","<<bx.yu<<"]"; }

struct PlanarProjectionMap { uint n, i, j; PlanarProjectionMap(uint nn, uint ii, uint jj) : n(nn), i(ii), j(jj) { } uint argument_size() const { return n; } };
inline std::ostream& operator<<(std::ostream& os, const PlanarProjectionMap& p) {
    return os << "PlanarProjectionMap( argument_size="<<p.n<<", x="<<p.i<<", y="<<p.j<<" )";
}

//! \brief Base interface for plotting and drawing classes.
class FigureInterface {
  public:
    virtual ~FigureInterface() { };
    virtual void set_projection_map(const PlanarProjectionMap& prj) = 0;
    virtual void set_bounding_box(const Box& bx) = 0;
    virtual void set_projection(uint as, uint ix, uint iy) = 0;
    virtual void set_line_style(bool) = 0;
    virtual void set_line_width(double) = 0;
    virtual void set_line_colour(Colour) = 0;
    virtual void set_fill_opacity(double) = 0;
    virtual void set_fill_colour(Colour) = 0;
    virtual bool get_line_style() const = 0;
    virtual double get_line_width() const = 0;
    virtual Colour get_line_colour() const = 0;
    virtual bool get_fill_style() const = 0;
    virtual double get_fill_opacity() const = 0;
    virtual Colour get_fill_colour() const = 0;
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
    virtual uint x_size_in_pixels() const = 0;
    virtual uint y_size_in_pixels() const = 0;
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
    virtual void set_fill_opacity(double fo) = 0;
    virtual void set_fill_colour(double r, double g, double b) = 0;
    virtual void get_bounding_box(double& xl, double& xu, double& yl, double& yu) const = 0;
  public:
    template<class X, class Y> void move_to(X x, Y y) { this->move_to(numeric_cast<double>(x),numeric_cast<double>(y)); }
    template<class X, class Y> void line_to(X x, Y y) { this->line_to(numeric_cast<double>(x),numeric_cast<double>(y)); }
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
