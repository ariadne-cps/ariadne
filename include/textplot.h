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
class TaylorSet;
class GridTreeSet;

    
//! \brief Class for plotting sets as a list of points.
class TextPlot 
    : public GraphicsInterface
{
  public:
    ~TextPlot();
    TextPlot();
    TextPlot(const char* filename);
    TextPlot(const char* filename, ios::openmode mode);
    void open(const char* filename);
    void open(const char* filename, ios::openmode mode);
    void draw(const std::vector<Point>&); // Draw a shape bounded by a list of points
    void draw(const Point&);
    void draw(const Box&);
    void draw(const Polytope&);
    void draw(const InterpolatedCurve&);
    void close();
  private:
    std::ofstream _fstream;
};

template<class SET> TextPlot& operator<<(TextPlot& g, const SET& set) { draw(g, set); return g; }

template<class SET> void textplot(const char* filename, const SET& set) { 
    TextPlot g(filename); draw(g, set); g.close(); }

} // namespace Ariadne

#endif // ARIADNE_TEXTPLOT_H
