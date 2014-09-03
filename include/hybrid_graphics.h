/***************************************************************************
 *            hybrid_graphics.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file hybrid_graphics.h
 *  \brief Graphics class for drawing objects in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_GRAPHICS_H
#define ARIADNE_HYBRID_GRAPHICS_H

#include <iosfwd>
#include <string>
#include <vector>

#include "colour.h"
#include "variables.h"
#include "graphics_interface.h"
#include "graphics.h"

#include "expression_set.h"

typedef unsigned int uint;

namespace Ariadne {

class IntervalSet;

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

Interval approximate_interval(const RealVariableInterval&);


struct Axes2d {
    Axes2d(const RealVariableInterval x, const RealVariableInterval& y)
            : variables(x.variable(),y.variable()), bounds() {
        bounds.insert(x.variable(),x.approximate_interval());
        bounds.insert(y.variable(),y.approximate_interval()); }
    Axes2d(double xl, const RealVariable& x, double xu, double yl, const RealVariable& y, double yu)
            : variables(x,y), bounds() {
        bounds.insert(x,Interval(xl,xu));
        bounds.insert(y,Interval(yl,yu)); }
    Variables2d variables;
    Map<RealVariable,Interval> bounds;
};

//! \brief Class for plotting figures of hybrid sets.
class HybridFigure
{
  public:
    ~HybridFigure();
    HybridFigure();

    void set_locations(const List<DiscreteLocation>& l) { locations=Set<DiscreteLocation>(l); }
    void set_axes(const Axes2d& axes) { bounds=axes.bounds; variables=axes.variables; }
    void set_bounds(const RealVariable& x, const Float& l, const Float& u) { bounds.insert(x,Interval(l,u)); }
    void set_bounds(const RealVariable& x, const Interval& ivl) { bounds.insert(x,ivl); }
    void set_bounds(const Map<RealVariable,Interval>& b) { bounds=b; };
    void set_bounds(const Map<RealVariable,IntervalSet>& b);
    void set_variables(const RealVariable& x, const RealVariable& y) { variables=Variables2d(x,y); }

    void set_line_style(bool ls) { properties.line_style=ls; }
    void set_line_width(double lw) { properties.line_width=lw; }
    void set_line_colour(Colour lc) { properties.line_colour=lc; }
    void set_fill_style(bool fs) { properties.fill_style=fs; }
    void set_fill_colour(Colour fc) { properties.fill_colour=fc; }

    void set_fill_opacity(double fo) { properties.fill_colour.opacity=fo; }
    void set_line_colour(double r, double g, double b) { properties.line_colour=Colour(r,g,b); }
    void set_fill_colour(double r, double g, double b) { properties.fill_colour=Colour(r,g,b,properties.fill_colour.opacity); }

    bool get_line_style() const { return properties.line_style; }
    double get_line_width() const { return properties.line_width; }
    Colour get_line_colour() const { return properties.line_colour; }
    bool get_fill_style() const { return properties.fill_style; }
    Colour get_fill_colour() const { return properties.fill_colour; }

    void draw(const HybridDrawableInterface& shape) { objects.append(HybridGraphicsObject(this->properties,shape)); }
    void clear() { objects.clear(); }

    void write(const char* filename, uint nx, uint ny) const;
    void write(const char* filename) const;
  public:
    void _paint_all(CanvasInterface& canvas) const; // Writes all shapes to the canvas
  private:
  public:
    Map<RealVariable,Interval> bounds;
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

inline void draw(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); }
inline HybridFigure& operator<<(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); return fig; }

Interval approximation(const IntervalSet& rivl);

template<class SET1>
void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1) {
    HybridFigure g; g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.write(filename); }

template<class SET1,class SET2>
void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2) {
    HybridFigure g; g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2); g.write(filename); }

template<class SET1,class SET2,class SET3>
void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3) {
    HybridFigure g; g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.write(filename); }

template<class SET1,class SET2,class SET3,class SET4>
void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4) {
    HybridFigure g;  g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4); g.write(filename); }

template<class SET1,class SET2,class SET3,class SET4,class SET5>
void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4, const Colour& fc5, const SET5& set5) {
    HybridFigure g;  g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4);
    g.set_fill_colour(fc5); draw(g,set5); g.write(filename); }

template<class SET1,class SET2,class SET3,class SET4,class SET5,class SET6>
void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4,
          const Colour& fc5, const SET5& set5, const Colour& fc6, const SET6& set6) {
    HybridFigure g;  g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4);
    g.set_fill_colour(fc5); draw(g,set5); g.set_fill_colour(fc6); draw(g,set6); g.write(filename); }

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
