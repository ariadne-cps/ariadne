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

#include "graphics_interface.h"

typedef unsigned int uint;

namespace Ariadne {

using namespace std;

class Point;
class Box;
class Polytope;
class InterpolatedCurve;
class Zonotope;
class TaylorImageSet;
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

    void set_projection(uint, uint, uint) { };

    void set_line_style(bool) { };
    void set_line_width(double) { };
    void set_line_colour(Colour) { };
    void set_fill_style(bool) { };
    void set_fill_opacity(double) { };
    void set_fill_colour(Colour) { };

    void set_line_colour(double, double, double) { };
    void set_fill_colour(double, double, double) { };

    bool get_line_style() const { return true; }
    double get_line_width() const { return 1.0; }
    Colour get_line_colour() const { return black; }
    bool get_fill_style() const { return false; };
    double get_fill_opacity() const { return 0.0; };
    Colour get_fill_colour() const { return white; };

    void open(const char* filename);
    void open(const char* filename, ios::openmode mode);
    void draw(const Point&);
    void draw(const Box&);
    void draw(const Polytope&);
    void draw(const InterpolatedCurve&);
    void draw(const GridTreeSubset&);
    void draw(const DrawableInterface&);
    void close();
  private:
    void _draw(const std::vector<Point>&);
  private:
    std::ofstream _fstream;
};

template<class SET> TextPlot& operator<<(TextPlot& g, const SET& set) { draw(g, set); return g; }

template<class SET> void textplot(const char* filename, const SET& set) {
    TextPlot g(filename); draw(g, set); g.close(); }

} // namespace Ariadne

#endif // ARIADNE_TEXTPLOT_H
