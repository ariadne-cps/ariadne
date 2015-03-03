/***************************************************************************
 *            textplot.h
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

/*! \file textplot.h
 *  \brief TextPlot class for outputting sets as a list of point (that can be imported in GnuPlot, Matlab, etc.).
 */

#ifndef ARIADNE_TEXTPLOT_H
#define ARIADNE_TEXTPLOT_H

#include <iosfwd>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "output/graphics_interface.h"
#include "output/colour.h"

typedef unsigned int Nat;

namespace Ariadne {

using namespace std;

template<class X> class Point;
typedef Point<ExactNumericType> ExactPoint;
template<class IVL> class Box;
typedef Box<ExactIntervalType> ExactBoxType;

class InterpolatedCurve;
class Zonotope;
class TaylorConstrainedImageSet;
class GridTreeSubset;


//! \brief Class for plotting sets as a list of points.
class TextPlot
    : public FigureInterface
{
  public:
    ~TextPlot();
    TextPlot();
    TextPlot(const char* filename);
    TextPlot(const char* filename, ios::openmode mode);

    virtual Void set_projection_map(const PlanarProjectionMap& prj) override { };
    virtual Void set_bounding_box(const GraphicsBoundingBoxType& bx) override { };
    virtual Void set_dot_radius(double) override { };

    Void set_projection(Nat, Nat, Nat) { };

    Void set_line_style(Bool) { };
    Void set_line_width(double) { };
    Void set_line_colour(Colour) { };
    Void set_fill_style(Bool) { };
    Void set_fill_opacity(double) { };
    Void set_fill_colour(Colour) { };

    Void set_line_colour(double, double, double) { };
    Void set_fill_colour(double, double, double) { };

    Bool get_line_style() const { return true; }
    double get_line_width() const { return 1.0; }
    Colour get_line_colour() const { return black; }
    Bool get_fill_style() const { return false; };
    double get_fill_opacity() const { return 0.0; };
    Colour get_fill_colour() const { return white; };

    Void open(const char* filename);
    Void open(const char* filename, ios::openmode mode);
    Void draw(const ExactPoint&);
    Void draw(const ExactBoxType&);
//    Void draw(const Polytope&);
    Void draw(const InterpolatedCurve&);
    Void draw(const GridTreeSubset&);
    Void draw(const DrawableInterface&);
    Void close();
  private:
    Void _draw(const std::vector<ExactPoint>&);
  private:
    std::ofstream _fstream;
};

template<class SET> TextPlot& operator<<(TextPlot& g, const SET& set) { draw(g, set); return g; }

template<class SET> Void textplot(const char* filename, const SET& set) {
    TextPlot g(filename); draw(g, set); g.close(); }

} // namespace Ariadne

#endif // ARIADNE_TEXTPLOT_H
