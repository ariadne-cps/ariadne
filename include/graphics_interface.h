/***************************************************************************
 *            graphics_interface.h
 *
 *  Copyright 2009-10  Davide Bresolin, Pieter Collins
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

typedef unsigned int uint;

namespace Ariadne {

template<class R, class A> inline R numeric_cast(const A&);
template<class T> class Set;

struct Vector2d;
struct Point2d;
struct Box2d;

class ApproximateInterval;
template<class IVL> class Box;
typedef Box<ApproximateInterval> ApproximateBox;

typedef ApproximateBox GraphicsBoundingBoxType;
class Colour;

class DrawableInterface;
class FigureInterface;
class CanvasInterface;

struct PlanarProjectionMap {
    uint n, i, j;
    PlanarProjectionMap(uint nn, uint ii, uint jj) : n(nn), i(ii), j(jj) { }
    uint argument_size() const { return n; }
    uint x_coordinate() const { return i; }
    uint y_coordinate() const { return j; }
};
inline std::ostream& operator<<(std::ostream& os, const PlanarProjectionMap& p) {
    return os << "P<R"<<p.n<<";R2>[x"<<p.i<<",x"<<p.j<<"]"; }
typedef PlanarProjectionMap Projection2d;

//! \ingroup GraphicsModule
//! \brief Base interface for plotting and drawing classes.
class FigureInterface {
  public:
    virtual ~FigureInterface() { };
    virtual void set_projection_map(const PlanarProjectionMap& prj) = 0;
    virtual void set_bounding_box(const ApproximateBox& bx) = 0;
    virtual void set_projection(uint as, uint ix, uint iy) = 0;
    virtual void set_line_style(bool) = 0;
    virtual void set_line_width(double) = 0;
    virtual void set_dot_radius(double) = 0;
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
inline void draw(FigureInterface& fig, const DrawableInterface& shape) { fig.draw(shape); }
inline FigureInterface& operator<<(FigureInterface& fig, const DrawableInterface& shape) { fig.draw(shape); return fig; }

class Figure;
void draw(Figure& fig, const DrawableInterface& shape);
Figure& operator<<(Figure& fig, const DrawableInterface& shape);

//! \ingroup GraphicsModule
//! \brief Interface to two-dimensional drawing canvas with the ability to draw polyhedra.
class CanvasInterface {
  public:
    //! \brief Destructor
    virtual ~CanvasInterface() { };

    virtual void initialise(std::string x, std::string y, double lx, double ux, double ly, double uy) = 0;
    virtual void finalise() = 0;

    //! \brief Move the current initial point for a line to the point \a (x,y).
    virtual void move_to(double x, double y) = 0;
    //! \brief Create a line segment from the current point to the point \a (x,y).
    virtual void line_to(double x, double y) = 0;
    //! \brief Draw a circle with centre \a (x,y) and radius \a r.
    virtual void circle(double x, double y, double r) = 0;
    //! \brief Draw a dot with centre \a (x,y).
    virtual void dot(double x, double y) = 0;
    //! \brief Draw the working line.
    virtual void stroke() = 0;
    //! \brief Draw and fill the working line.
    virtual void fill() = 0;
    //! Set the width of the lines bounding shapes, in pixels.
    virtual void set_line_width(double lw) = 0;
    //! \brief Set the colour of subsequent line to be drawn. The \a r, \a g and \a b
    //! arguments are intensities of red, green and blue on a scale of 0.0 to 1.0.
    virtual void set_line_colour(double r, double g, double b) = 0;
    //! \brief  Set the colour of the next regions to be filled.
    virtual void set_fill_opacity(double fo) = 0;
    //! \brief Set the colour of subsequent regions to be filled.
    virtual void set_fill_colour(double r, double g, double b) = 0;

    //! \brief The scaling of the figure, in user units per pixel.
    virtual Vector2d scaling() const = 0;
    //! brief The lower and upper bounds of the x- and y- coordinates of the drawing region.
    virtual Box2d bounds() const = 0;
  public:
    template<class X, class Y> void move_to(X x, Y y) { this->move_to(numeric_cast<double>(x),numeric_cast<double>(y)); }
    template<class X, class Y> void line_to(X x, Y y) { this->line_to(numeric_cast<double>(x),numeric_cast<double>(y)); }
};


//! \ingroup GraphicsModule
//! \brief Base interface for drawable objects
class DrawableInterface {
  public:
    //! brief Virtual destructor.
    virtual ~DrawableInterface() { }
    //! brief Make a dynamically-allocated copy.
    virtual DrawableInterface* clone() const = 0;
    //! brief Draw the object on the canvas \a c using line segments and fill/stroke commands.
    virtual void draw(CanvasInterface& c, const Projection2d& p) const = 0;
    //! brief The dimension of the object in Euclidean space
    virtual uint dimension() const = 0;
    //! brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const { return os << "Drawable"; }
};




class Real;
template<class T> class Variable;
typedef Variable<Real> RealVariable;
template<class T> class Space;
typedef Space<Real> RealSpace;
class DiscreteLocation;

struct Variables2d;

bool valid_axis_variables(const RealSpace& space, const Variables2d& variables);
Projection2d projection(const RealSpace& spc, const Variables2d& variables);

//! \ingroup GraphicsModule
//! \brief Base interface for drawable objects
class HybridDrawableInterface {
  public:
    //! brief Virtual destructor.
    virtual ~HybridDrawableInterface() { }
    //! brief Draw the object on the canvas \a c using line segments and fill/stroke commands.
    virtual void draw(CanvasInterface& c, const Set<DiscreteLocation>& q, const Variables2d& v) const = 0;
};


} // namespace Ariadne


#endif // ARIADNE_GRAPHICS_INTERFACE_H
