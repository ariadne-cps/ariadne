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

class RealInterval;

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

bool valid_axes(const RealSpace& space, const Variables2d& axes);
Projection2d projection(const RealSpace& spc, const Variables2d& axes);

//! \brief Class for plotting figures of hybrid sets.
class HybridFigure
{
  public:
    ~HybridFigure();
    HybridFigure();

    void set_locations(const List<DiscreteLocation>& l) { locations=Set<DiscreteLocation>(l); }
    void set_bounds(const RealVariable& x, const Float& l, const Float& u) { bounds.insert(x,Interval(l,u)); }
    void set_bounds(const RealVariable& x, const Interval& ivl) { bounds.insert(x,ivl); }
    void set_bounds(const Map<RealVariable,Interval>& b) { bounds=b; };
    void set_bounds(const Map<RealVariable,RealInterval>& b);
    void set_axes(const RealVariable& x, const RealVariable& y) { axes=Variables2d(x,y); }

    void set_line_style(bool ls) { properties.line_style=ls; }
    void set_line_width(double lw) { properties.line_width=lw; }
    void set_line_colour(Colour lc) { properties.line_colour=lc; }
    void set_fill_style(bool fs) { properties.fill_style=fs; }
    void set_fill_opacity(double fo) { properties.fill_opacity=fo; }
    void set_fill_colour(Colour fc) { properties.fill_colour=fc; }

    void set_line_colour(double r, double g, double b) { properties.line_colour=Colour(r,g,b); }
    void set_fill_colour(double r, double g, double b) { properties.fill_colour=Colour(r,g,b); }

    bool get_line_style() const { return properties.line_style; }
    double get_line_width() const { return properties.line_width; }
    Colour get_line_colour() const { return properties.line_colour; }
    bool get_fill_style() const { return properties.fill_style; }
    double get_fill_opacity() const { return properties.fill_opacity; }
    Colour get_fill_colour() const { return properties.fill_colour; }

    void draw(const HybridDrawableInterface& shape) { objects.append(HybridGraphicsObject(this->properties,shape)); }
    void clear() { objects.clear(); }

    void write(const char* filename, uint nx, uint ny) const;
    void write(const char* filename) const;
  public:
    void _paint_all(CanvasInterface& canvas) const; // Writes all shapes to the canvas
  private:
    Map<RealVariable,Interval> bounds;
    Set<DiscreteLocation> locations;
    Variables2d axes;
    GraphicsProperties properties;
    List<HybridGraphicsObject> objects;
};

inline HybridFigure& operator<<(HybridFigure& g, const LineStyle& ls) { g.set_line_style(ls); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const LineWidth& lw) { g.set_line_width(lw); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const LineColour& lc) { g.set_line_colour(lc); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const FillStyle& fs) { g.set_fill_style(fs); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const FillOpacity& fo) { g.set_fill_opacity(fo); return g; }
inline HybridFigure& operator<<(HybridFigure& g, const FillColour& fc) { g.set_fill_colour(fc); return g; }

inline void hdraw(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); }
inline HybridFigure& operator<<(HybridFigure& fig, const HybridDrawableInterface& shape) { fig.draw(shape); return fig; }

Interval approximation(const RealInterval& rivl);

template<class SET1>
void hplot(const char* filename, const List<RealVariableInterval>& bbox, const Colour& fc1, const SET1& set1) {
    HybridFigure g; for(uint i=0; i!=bbox.size(); ++i) { g.set_bounds(bbox[i].variable(),approximation(bbox[i].interval())); }
    g.set_axes(bbox[0].variable(),bbox[1].variable()); g.set_fill_colour(fc1); hdraw(g,set1); g.write(filename); }

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
