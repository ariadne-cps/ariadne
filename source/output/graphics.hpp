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
    friend OutputStream& operator<<(OutputStream& os, GraphicsProperties const& gp);
};

Void set_properties(CanvasInterface& canvas, const GraphicsProperties& properties);

//! \brief Class for plotting figures.
class Figure
    : public FigureInterface
{
  public:
    ~Figure();
    Figure(); //< Deprecated
    Figure(const GraphicsBoundingBoxType& bbx, const PlanarProjectionMap& proj);
    //! Construct a figure projecting \a bbx onto the (\a ix, \a iy) coordinates
    Figure(const GraphicsBoundingBoxType& bbx, DimensionType ix, DimensionType iy);
    //! Construct a figure projecting \a bbx onto the (\a ix, \a iy) coordinates
    Figure(const GraphicsBoundingBoxType& bbx, Pair<DimensionType,DimensionType> ixy);
    Figure& set_projection_map(const PlanarProjectionMap&);
    //! Set the restricted region to display coordinates to (\a i, \a j)
    Figure& set_bounding_box(const GraphicsBoundingBoxType&);

    PlanarProjectionMap get_projection_map() const;
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

inline Figure& operator<<(Figure& g, const LineStyle& ls) { g.set_line_style(ls); return g; }
inline Figure& operator<<(Figure& g, const LineWidth& lw) { g.set_line_width(lw); return g; }
inline Figure& operator<<(Figure& g, const LineColour& lc) { g.set_line_colour(lc); return g; }
inline Figure& operator<<(Figure& g, const FillStyle& fs) { g.set_fill_style(fs); return g; }
inline Figure& operator<<(Figure& g, const FillOpacity& fo) { g.set_fill_opacity(fo); return g; }
inline Figure& operator<<(Figure& g, const FillColour& fc) { g.set_fill_colour(fc); return g; }

inline Figure& operator<<(Figure& fig, const DrawableInterface& shape) { fig.draw(shape); return fig; }
inline Figure& operator<<(Figure& fig, const ApproximateBoxType& box) { fig.draw(box); return fig; }

inline Void draw(Figure& g) { }

template<class SET, class... CSETS> inline Void
draw(Figure& g, const Colour& fc, const SET& set, CSETS const& ... csets) {
    g.set_fill_colour(fc); draw(g,set); draw(g,csets...); }

template<class... CSETS> Void
plot(const char* filename, const PlanarProjectionMap& pr, const ApproximateBoxType& bbox, CSETS const&... csets) {
    Figure g; g.set_projection_map(pr); g.set_bounding_box(bbox); draw(g,csets...);  g.write(filename); }

template<class... CSETS> Void
plot(const char* filename, const ApproximateBoxType& bbox, CSETS const&... csets) {
    plot(filename, PlanarProjectionMap(2u,0,1), bbox, csets...); }

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_HPP
