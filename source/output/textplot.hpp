/***************************************************************************
 *            textplot.hpp
 *
 *  Copyright 2009-17  Davide Bresolin
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

/*! \file textplot.hpp
 *  \brief TextPlot class for outputting sets as a list of point (that can be imported in GnuPlot, Matlab, etc.).
 */

#ifndef ARIADNE_TEXTPLOT_HPP
#define ARIADNE_TEXTPLOT_HPP

#include <iosfwd>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "../output/graphics_interface.hpp"
#include "../output/colour.hpp"

namespace Ariadne {

typedef unsigned int Nat;

template<class X> class Point;
typedef Point<ExactNumericType> ExactPoint;
template<class IVL> class Box;
typedef Box<ExactIntervalType> ExactBoxType;

class InterpolatedCurve;
class Zonotope;
class TaylorConstrainedImageSet;
class GridTreeSubpaving;


//! \brief Class for plotting sets as a list of points.
class TextPlot
    : public FigureInterface
{
  public:
    ~TextPlot();
    TextPlot();
    TextPlot(const char* filename);
    TextPlot(const char* filename, std::ios::openmode mode);

    virtual TextPlot& set_projection_map(const PlanarProjectionMap& prj) override { return *this; };
    virtual TextPlot& set_bounding_box(const GraphicsBoundingBoxType& bx) override { return *this; };
    virtual TextPlot& set_dot_radius(double) override { return *this; };

    TextPlot& set_projection(Nat, Nat, Nat) override { return *this; };

    TextPlot& set_line_style(Bool) override { return *this; };
    TextPlot& set_line_width(double) override { return *this; };
    TextPlot& set_line_colour(Colour) override { return *this; };
    TextPlot& set_fill_opacity(double) override { return *this; };
    TextPlot& set_fill_colour(Colour) override { return *this; };

    TextPlot& set_line_colour(double, double, double) { return *this; };
    TextPlot& set_fill_colour(double, double, double) { return *this; };

    Bool get_line_style() const override { return true; }
    double get_line_width() const override { return 1.0; }
    Colour get_line_colour() const override { return black; }
    Bool get_fill_style() const override { return false; };
    double get_fill_opacity() const override { return 0.0; };
    Colour get_fill_colour() const override { return white; };

    Void open(const char* filename);
    Void open(const char* filename, std::ios::openmode mode);
    TextPlot& draw(const ExactPoint&);
    TextPlot& draw(const ExactBoxType&);
//    TextPlot& draw(const Polytope&);
    TextPlot& draw(const InterpolatedCurve&);
    TextPlot& draw(const GridTreeSubpaving&);
    TextPlot& draw(const DrawableInterface&) override;
    Void close();
  private:
    TextPlot& _draw(const std::vector<ExactPoint>&);
  private:
    std::ofstream _fstream;
};

template<class SET> TextPlot& operator<<(TextPlot& g, const SET& set) { draw(g, set); return g; }

template<class SET> Void textplot(const char* filename, const SET& set) {
    TextPlot g(filename); draw(g, set); g.close(); }

} // namespace Ariadne

#endif // ARIADNE_TEXTPLOT_HPP
