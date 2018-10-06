/***************************************************************************
 *            hybrid_graphics.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file hybrid_graphics.hpp
 *  \brief Graphics class for drawing objects in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_GRAPHICS_HPP
#define ARIADNE_HYBRID_GRAPHICS_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "../output/colour.hpp"
#include "../symbolic/variables.hpp"
#include "../output/graphics_interface.hpp"
#include "../output/graphics.hpp"

#include "../symbolic/expression_set.hpp"
#include "../hybrid/discrete_location.hpp"
#include "../hybrid/hybrid_graphics_interface.hpp"

namespace Ariadne {

template<> class Interval<FloatDP> {
    FloatDP _l, _u;
  public:
    Interval(FloatDP l, FloatDP u) : _l(l), _u(u) { }
    template<class UB> Interval(Interval<UB> const& ivl) : _l(ivl.lower()), _u(ivl.upper()) { }
    FloatDP lower() const { return _l; }
    FloatDP upper() const { return _u; }
};

template<> class Interval<ApproximateDouble> {
    ApproximateDouble _l, _u;
  public:
    Interval(ApproximateDouble l, ApproximateDouble u) : _l(l), _u(u) { }
    template<class UB> Interval(Interval<UB> const& ivl) : _l(ivl.lower()), _u(ivl.upper()) { }
    ApproximateDouble lower() const { return _l; }
    ApproximateDouble upper() const { return _u; }
};
using ApproximateDoubleInterval = Interval<ApproximateDouble>;
using ApproximateDoubleVariableInterval = VariableInterval<ApproximateDouble>;


struct HybridGraphicsObject {
    HybridGraphicsObject(const GraphicsProperties& gp, const HybridDrawableInterface& sh)
        : properties(gp), shape_ptr(&sh) { }
    GraphicsProperties properties;
    const HybridDrawableInterface* shape_ptr;
};

struct Variables2d {
    Identifier _x,_y;
    Variables2d(const RealVariable& x, const RealVariable& y) : _x(x.name()), _y(y.name()) { }
    RealVariable x_variable() const { return RealVariable(_x); }
    RealVariable y_variable() const { return RealVariable(_y); }
};



struct Axes2d {
    Axes2d(const ApproximateDoubleVariableInterval x, const ApproximateDoubleVariableInterval& y);
    Axes2d(ApproximateDouble xl, const RealVariable& x, ApproximateDouble xu, ApproximateDouble yl, const RealVariable& y, ApproximateDouble yu);
    Variables2d variables;
    Map<RealVariable,ApproximateDoubleInterval> bounds;
};

inline Axes2d::Axes2d(const ApproximateDoubleVariableInterval x, const ApproximateDoubleVariableInterval& y)
        : variables(x.variable(),y.variable()), bounds() {
    bounds.insert(x.variable(),x.interval());
    bounds.insert(y.variable(),y.interval());
}

inline Axes2d::Axes2d(ApproximateDouble xl, const RealVariable& x, ApproximateDouble xu, ApproximateDouble yl, const RealVariable& y, ApproximateDouble yu)
        : variables(x,y), bounds() {
    bounds.insert(x,ApproximateDoubleInterval(xl,xu));
    bounds.insert(y,ApproximateDoubleInterval(yl,yu));
}

//! \brief Class for plotting figures of hybrid sets.
class HybridFigure
{
  public:
    ~HybridFigure();
    HybridFigure();

    Void set_locations(const List<DiscreteLocation>& l) { locations=Set<DiscreteLocation>(l); }
    Void set_axes(const Axes2d& axes) { bounds=axes.bounds; variables=axes.variables; }
    Void set_bounds(const RealVariable& x, const ApproximateDouble& l, const ApproximateDouble& u) { bounds.insert(x,ApproximateDoubleInterval(l,u)); }
    Void set_bounds(const RealVariable& x, const ApproximateDoubleInterval& ivl) { bounds.insert(x,ivl); }
    Void set_bounds(const Map<RealVariable,ApproximateDoubleInterval>& b) { bounds=b; };
    Void set_variables(const RealVariable& x, const RealVariable& y) { variables=Variables2d(x,y); }

    Void set_line_style(Bool ls) { properties.line_style=ls; }
    Void set_line_width(double lw) { properties.line_width=lw; }
    Void set_line_colour(Colour lc) { properties.line_colour=lc; }
    Void set_fill_style(Bool fs) { properties.fill_style=fs; }
    Void set_fill_colour(Colour fc) { properties.fill_colour=fc; }

    Void set_fill_opacity(double fo) { properties.fill_colour.opacity=fo; }
    Void set_line_colour(double r, double g, double b) { properties.line_colour=Colour(r,g,b); }
    Void set_fill_colour(double r, double g, double b) { properties.fill_colour=Colour(r,g,b,properties.fill_colour.opacity); }

    Bool get_line_style() const { return properties.line_style; }
    double get_line_width() const { return properties.line_width; }
    Colour get_line_colour() const { return properties.line_colour; }
    Bool get_fill_style() const { return properties.fill_style; }
    Colour get_fill_colour() const { return properties.fill_colour; }

    Void draw(const HybridDrawableInterface& shape) { objects.append(HybridGraphicsObject(this->properties,shape)); }
    Void clear() { objects.clear(); }

    Void write(const char* filename, Nat nx, Nat ny) const;
    Void write(const char* filename) const;
  public:
    Void _paint_all(CanvasInterface& canvas) const; // Writes all shapes to the canvas
  private:
  public:
    Map<RealVariable,ApproximateDoubleInterval> bounds;
    Set<DiscreteLocation> locations;
    Variables2d variables;
    GraphicsProperties properties;
    List<HybridGraphicsObject> objects;
};

inline HybridFigure& operator<<(HybridFigure& g, const LineStyle& ls) { g.set_line_style(ls); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const LineWidth& lw) { g.set_line_width(lw); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const LineColour& lc) { g.set_line_colour(lc); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const FillStyle& fs) { g.set_fill_style(fs); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const FillOpacity& fo) { g.set_fill_opacity(fo); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const FillColour& fc) { g.set_fill_colour(fc); return g; }

inline Void draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& variables, const HybridDrawableInterface& shape) { shape.draw(canvas,locations,variables); }
inline Void draw(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); }

Void paint(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& variables, const List<HybridGraphicsObject>& objects);

inline HybridFigure& operator<<(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); return fig; }

inline Void draw(HybridFigure& g) { }

template<class SET, class... CSETS>
inline Void draw(HybridFigure& g, const Colour& fc1, const SET& set1, CSETS const&... csets) {
    g.set_fill_colour(fc1); draw(g,set1); draw(g,csets...); }

template<class... CSETS>
Void plot(const char* filename, const Axes2d& axes, CSETS const&... csets) {
    HybridFigure g;  g.set_axes(axes); draw(g,csets...); g.write(filename); }

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_HPP
