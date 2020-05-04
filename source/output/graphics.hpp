/***************************************************************************
 *            output/graphics.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file output/graphics.hpp
 *  \brief Graphics class for drawing and outputting shapes in Euclidean space.
 */

#ifndef ARIADNE_GRAPHICS_HPP
#define ARIADNE_GRAPHICS_HPP

#include "../config.hpp"

#include <iosfwd>
#include <string>
#include <vector>

#include "../utility/typedefs.hpp"
#include "../utility/declarations.hpp"
#include "../symbolic/variables.hpp"
#include "../output/colour.hpp"
#include "../output/graphics_interface.hpp"

namespace Ariadne {

typedef ApproximateBoxType GraphicsBoundingBoxType;

class ProjectionFunction;

struct LineStyle { explicit LineStyle(Bool ls) : _style(ls) { } operator Bool() const { return this->_style; } private: Bool _style; };
struct LineWidth { explicit LineWidth(Dbl lw) : _width(lw) { } operator Dbl() const { return this->_width; } private: Dbl _width; };
struct LineColour : Colour { LineColour(const Colour& lc) : Colour(lc) { } LineColour(Dbl r, Dbl g, Dbl b) : Colour(r,g,b) { } };
struct FillStyle { explicit FillStyle(Bool fs) : _style(fs) { } operator Bool() const { return this->_style; } private: Bool _style; };
struct FillOpacity { explicit FillOpacity(Dbl fo) : _opacity(fo) { } operator Dbl() const { return this->_opacity; } private: Dbl _opacity; };
struct FillColour : Colour { FillColour(const Colour& fc) : Colour(fc) { } FillColour(Dbl r, Dbl g, Dbl b) : Colour(r,g,b) { } };

inline LineStyle line_style(Bool s) { return LineStyle(s); }
inline LineWidth line_width(Dbl w) { return LineWidth(w); }
inline LineColour line_colour(const Colour& c) { return LineColour(c); }
inline LineColour line_colour(Dbl r, Dbl g, Dbl b) { return LineColour(Colour(r,g,b)); }
inline FillStyle fill_style(Bool s) { return FillStyle(s); }
inline FillOpacity fill_opacity(Dbl o) { return FillOpacity(o); }
inline FillColour fill_colour(const Colour& c) { return FillColour(c); }
inline FillColour fill_colour(Dbl r, Dbl g, Dbl b) { return FillColour(Colour(r,g,b)); }

struct GraphicsProperties {
    GraphicsProperties()
        : dot_radius(1.0), line_style(true), line_width(1.0), line_colour(black), fill_style(true), fill_colour(ariadneorange) { }
    GraphicsProperties(Bool ls, Dbl lw, Dbl dr, Colour lc, Bool fs, Colour fc)
        : dot_radius(dr), line_style(ls), line_width(lw), line_colour(lc), fill_style(fs), fill_colour(fc) { }
    Dbl dot_radius;
    Bool line_style;
    Dbl line_width;
    Colour line_colour;
    Bool fill_style;
    Colour fill_colour;

    GraphicsProperties& set_dot_radius(Dbl);
    GraphicsProperties& set_line_style(Bool);
    GraphicsProperties& set_line_width(Dbl);
    GraphicsProperties& set_line_colour(Colour);
    GraphicsProperties& set_fill_style(Bool);
    GraphicsProperties& set_fill_colour(Colour);

    GraphicsProperties& set_line_colour(Dbl, Dbl, Dbl);
    GraphicsProperties& set_fill_colour(Dbl, Dbl, Dbl);
    GraphicsProperties& set_fill_opacity(Dbl);

    friend OutputStream& operator<<(OutputStream& os, GraphicsProperties const& gp);
};

Void set_properties(CanvasInterface& canvas, const GraphicsProperties& properties);

// TODO: Move to interval.hpp as templated IntervalData classes
template<> class Interval<Double> {
    Double _l, _u;
  public:
    Interval(Double l, Double u) : _l(l), _u(u) { }
    Double lower() const { return _l; }
    Double upper() const { return _u; }
};

template<> class Interval<ApproximateDouble> {
    ApproximateDouble _l, _u;
  public:
    typedef ApproximateDouble UpperBoundType;
    Interval(ApproximateDouble l, ApproximateDouble u) : _l(l), _u(u) { }
    template<class UB> Interval(Interval<UB> const& ivl) : _l(ivl.lower()), _u(ivl.upper()) { }
    ApproximateDouble lower() const { return _l; }
    ApproximateDouble upper() const { return _u; }
};
using ApproximateDoubleInterval = Interval<ApproximateDouble>;
using ApproximateDoubleVariableInterval = VariableInterval<ApproximateDouble>;
using DoubleVariableInterval = VariableInterval<Double>;

struct Variables2d {
    Identifier _x,_y;
    Variables2d(const RealVariable& x, const RealVariable& y);
    RealVariable x() const;
    RealVariable y() const;
    RealVariable x_variable() const;
    RealVariable y_variable() const;
};

struct Axes2d {
    Axes2d(const ApproximateDoubleVariableInterval x, const ApproximateDoubleVariableInterval& y);
    Axes2d(ApproximateDouble xl, const RealVariable& x, ApproximateDouble xu, ApproximateDouble yl, const RealVariable& y, ApproximateDouble yu);
    Variables2d variables;
    Map<RealVariable,ApproximateDoubleInterval> bounds;
};



//! \brief Class for plotting figures.
class Figure
    : public FigureInterface
{
  public:
    ~Figure();
    Figure(); //< Deprecated
    //! Construct a figure projecting \a bbx onto the \a proj coordinates
    Figure(const GraphicsBoundingBoxType& bbx, const Projection2d& proj);
    //! Construct a figure projecting \a bbx onto the (\a ix, \a iy) coordinates
    Figure(const GraphicsBoundingBoxType& bbx, DimensionType ix, DimensionType iy);
    //! Construct a figure projecting \a bbx onto the (\a ix, \a iy) coordinates
    Figure(const GraphicsBoundingBoxType& bbx, Pair<DimensionType,DimensionType> ixy);

    Figure& set_projection_map(const Projection2d&);
    Projection2d get_projection_map() const;
    Figure& set_bounding_box(const GraphicsBoundingBoxType&);
    GraphicsBoundingBoxType get_bounding_box() const;
    //! Set the displayed coordinates to (\a i, \a j)
    Figure& set_projection(DimensionType i, DimensionType j);
    Figure& set_projection(DimensionType as, DimensionType ix, DimensionType iy);

    //! Set the radiues to draw points (dots).
    Figure& set_dot_radius(Dbl);
    Figure& set_line_style(Bool);
    //! Set the width to draw lines. A width of 0 means no lines are drawn.
    Figure& set_line_width(Dbl);
    //! Set the colour draw lines
    Figure& set_line_colour(Colour);
    Figure& set_fill_style(Bool);
    //! Set the colour to fill shapes.
    Figure& set_fill_colour(Colour);

    Figure& set_line_colour(Dbl, Dbl, Dbl);
    Figure& set_fill_colour(Dbl, Dbl, Dbl);
    //! Set the opacity of shapes. An opacity of 0 means no fill.
    Figure& set_fill_opacity(Dbl);

    Bool get_line_style() const;
    Dbl get_line_width() const;
    Dbl get_dot_radius() const;
    Colour get_line_colour() const;
    Bool get_fill_style() const;
    Dbl get_fill_opacity() const;
    Colour get_fill_colour() const;

    GraphicsProperties& properties();
    GraphicsProperties const& properties() const;

    //! Add a set to draw onto the figure.
    Figure& draw(const DrawableInterface& shape);
    Figure& draw(const ApproximateBoxType& box);
    Figure& draw(const RealBox& box);

    //! Clear the figure.
    Figure& clear();
    //! Display the figure.
    Void display() const;
    //! Write out to file, using width \a nx pixels, and height \a ny pixels
    Void write(const Char* filename, Nat nx, Nat ny) const;
    //! Write to \a filename.
    Void write(const Char* filename) const;
  public:
    struct Data;
  public:
    Void _paint_all(CanvasInterface& canvas) const; // Writes all shapes to the canvas
  private:
    Data* _data;
};


Void draw(Figure& fig, const ApproximateBoxType& box);
inline Figure& operator<<(Figure& fig, const ApproximateBoxType& box) { fig.draw(box); return fig; }

inline Figure& operator<<(Figure& g, const LineStyle& ls) { g.properties().set_line_style(ls); return g; }
inline Figure& operator<<(Figure& g, const LineWidth& lw) { g.properties().set_line_width(lw); return g; }
inline Figure& operator<<(Figure& g, const LineColour& lc) { g.properties().set_line_colour(lc); return g; }
inline Figure& operator<<(Figure& g, const FillStyle& fs) { g.properties().set_fill_style(fs); return g; }
inline Figure& operator<<(Figure& g, const FillOpacity& fo) { g.properties().set_fill_opacity(fo); return g; }
inline Figure& operator<<(Figure& g, const FillColour& fc) { g.properties().set_fill_colour(fc); return g; }

template<class S> class LabelledSet;

//! \brief Class for plotting figures.
class LabelledFigure {
  public:
    ~LabelledFigure();

    //! Construct a figure drawing the given coordinates in the given bounds.
    LabelledFigure(const Axes2d& axes);

    LabelledFigure(const Variables2d& vars, VariablesBox<ApproximateIntervalType> const& bnds);

    Void set_axes(const Axes2d& axes);
//    Void set_bounds(const RealVariable& x, const ApproximateDouble& l, const ApproximateDouble& u);
    Void set_bounds(const RealVariable& x, const ApproximateDoubleInterval& ivl);
    Void set_bounds(const Map<RealVariable,ApproximateDoubleInterval>& b);

    Void set_bounding_box(VariablesBox<ApproximateIntervalType> const& bnds);

    GraphicsProperties& properties();
    GraphicsProperties const& properties() const;

    //! Add a set to draw onto the figure.
    LabelledFigure& draw(const LabelledDrawableInterface& shape);

    //! Clear the figure.
    LabelledFigure& clear();
    //! Display the figure.
    Void display() const;
    //! Write out to file, using width \a nx pixels, and height \a ny pixels
    Void write(const Char* filename, Nat nx, Nat ny) const;
    //! Write to \a filename.
    Void write(const Char* filename) const;
  public:
    struct Data;
  public:
    Void _paint_all(CanvasInterface& canvas) const; // Writes all shapes to the canvas
  private:
    Data* _data;
};

inline LabelledFigure& operator<<(LabelledFigure& g, const LineStyle& ls) { g.properties().set_line_style(ls); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const LineWidth& lw) { g.properties().set_line_width(lw); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const LineColour& lc) { g.properties().set_line_colour(lc); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const FillStyle& fs) { g.properties().set_fill_style(fs); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const FillOpacity& fo) { g.properties().set_fill_opacity(fo); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const FillColour& fc) { g.properties().set_fill_colour(fc); return g; }

template<class S> class LabelledSet;
template<class S> class LabelledDrawableWrapper : public LabelledDrawableInterface {
    LabelledSet<S> _lset;
  public:
    LabelledDrawableWrapper(LabelledSet<S> lset) : _lset(lset) { }
    virtual LabelledDrawableWrapper<S>* clone() const { return new LabelledDrawableWrapper<S>(*this); }
    virtual Void draw(CanvasInterface& cnvs, Variables2d const& vars) const {
        Projection2d prj(this->_lset.euclidean_set().dimension(),this->_lset.space()[vars.x()],this->_lset.space()[vars.y()]);
        this->_lset.euclidean_set().draw(cnvs,prj);
    }
};
template<class S> inline LabelledFigure& operator<<(LabelledFigure& g, const LabelledSet<S>& lset) {
    return g << LabelledDrawableWrapper<S>(lset);
}

inline Void draw(Figure& g) { }

template<class SET, class... CSETS> inline Void
draw(Figure& g, const Colour& fc, const SET& set, CSETS const& ... csets) {
    g.set_fill_colour(fc); draw(g,set); draw(g,csets...); }

template<class... CSETS> Void
plot(const char* filename, const Projection2d& pr, const ApproximateBoxType& bbox, CSETS const&... csets) {
    Figure g; g.set_projection_map(pr); g.set_bounding_box(bbox); draw(g,csets...);  g.write(filename); }

template<class... CSETS> Void
plot(const char* filename, const ApproximateBoxType& bbox, CSETS const&... csets) {
    plot(filename, Projection2d(2u,0,1), bbox, csets...); }

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_HPP

