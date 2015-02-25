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

#include "output/colour.h"
#include "expression/variables.h"
#include "output/graphics_interface.h"
#include "output/graphics.h"

#include "expression/expression_set.h"
#include "hybrid/discrete_location.h"

namespace Ariadne {

template<> class Interval<RawFloat64> {
    RawFloat64 _l, _u;
  public:
    Interval(RawFloat64 l, RawFloat64 u) :_l(l), _u(u) { }
    template<class UB> Interval(Interval<UB> const& ivl) : _l(ivl.lower()), _u(ivl.upper()) { }
    RawFloat64 lower() const { return _l; }
    RawFloat64 upper() const { return _u; }
};
typedef Interval<RawFloat64> FloatInterval;

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


struct FloatVariableLowerInterval {
    RawFloat64 _lower; RealVariable _variable;
    FloatVariableLowerInterval(const RawFloat64& l, const RealVariable& v) : _lower(l), _variable(v) { }
};

class FloatVariableInterval {
  private:
    RawFloat64 _lower; Variable<Real> _variable; RawFloat64 _upper;
  public:
    FloatVariableInterval(const RawFloat64& l, const Variable<Real>& v, const RawFloat64& u)
        : _lower(l), _variable(v), _upper(u) { ARIADNE_ASSERT_MSG(l<=u,"ExactInterval("<<l<<","<<u<<") not provably nonempty"); }
    FloatVariableInterval(const RealVariableInterval& rvivl)
        : _lower(rvivl.lower().get_d()), _variable(rvivl.variable()), _upper(rvivl.upper().get_d()) { }
    Variable<Real> const& variable() const { return this->_variable; }
    const FloatInterval interval() const { return FloatInterval(this->_lower,this->_upper); }
    const RawFloat64 lower() const { return this->_lower; }
    const RawFloat64 upper() const { return this->_upper; }
};
inline FloatVariableLowerInterval operator<=(double l, RealVariable const& v) {
    return FloatVariableLowerInterval(l,v); }
inline FloatVariableInterval operator<=(FloatVariableLowerInterval lv, double u) {
    return FloatVariableInterval(lv._lower,lv._variable,u); }
inline FloatVariableInterval operator<=(FloatVariableLowerInterval lv, RawFloat64 u) {
    return FloatVariableInterval(lv._lower,lv._variable,u); }
inline FloatVariableInterval operator<=(FloatVariableLowerInterval lv, Real u) {
    return FloatVariableInterval(lv._lower,lv._variable,RawFloat64(u.get_d())); }

struct Axes2d {
    Axes2d(const FloatVariableInterval x, const FloatVariableInterval& y)
            : variables(x.variable(),y.variable()), bounds() {
        bounds.insert(x.variable(),x.interval());
        bounds.insert(y.variable(),y.interval()); }
    Axes2d(double xl, const RealVariable& x, double xu, double yl, const RealVariable& y, double yu)
            : variables(x,y), bounds() {
        bounds.insert(x,FloatInterval(xl,xu));
        bounds.insert(y,FloatInterval(yl,yu)); }
    Variables2d variables;
    Map<RealVariable,FloatInterval> bounds;
};

//! \brief Class for plotting figures of hybrid sets.
class HybridFigure
{
  public:
    ~HybridFigure();
    HybridFigure();

    Void set_locations(const List<DiscreteLocation>& l) { locations=Set<DiscreteLocation>(l); }
    Void set_axes(const Axes2d& axes) { bounds=axes.bounds; variables=axes.variables; }
    Void set_bounds(const RealVariable& x, const RawFloat64& l, const RawFloat64& u) { bounds.insert(x,ExactInterval(l,u)); }
    Void set_bounds(const RealVariable& x, const FloatInterval& ivl) { bounds.insert(x,ivl); }
    Void set_bounds(const Map<RealVariable,FloatInterval>& b) { bounds=b; };
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
    Map<RealVariable,FloatInterval> bounds;
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

inline Void draw(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); }
inline HybridFigure& operator<<(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); return fig; }

template<class SET1>
Void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1) {
    HybridFigure g; g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.write(filename); }

template<class SET1,class SET2>
Void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2) {
    HybridFigure g; g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2); g.write(filename); }

template<class SET1,class SET2,class SET3>
Void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3) {
    HybridFigure g; g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.write(filename); }

template<class SET1,class SET2,class SET3,class SET4>
Void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4) {
    HybridFigure g;  g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4); g.write(filename); }

template<class SET1,class SET2,class SET3,class SET4,class SET5>
Void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4, const Colour& fc5, const SET5& set5) {
    HybridFigure g;  g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4);
    g.set_fill_colour(fc5); draw(g,set5); g.write(filename); }

template<class SET1,class SET2,class SET3,class SET4,class SET5,class SET6>
Void plot(const char* filename, const Axes2d& axes, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2,
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4,
          const Colour& fc5, const SET5& set5, const Colour& fc6, const SET6& set6) {
    HybridFigure g;  g.set_axes(axes); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2);
    g.set_fill_colour(fc3); draw(g,set3); g.set_fill_colour(fc4); draw(g,set4);
    g.set_fill_colour(fc5); draw(g,set5); g.set_fill_colour(fc6); draw(g,set6); g.write(filename); }

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
