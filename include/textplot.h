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
#include <map>

#include "graphics_interface.h"
#include "list_set.h"

typedef unsigned int uint;

namespace Ariadne {

using namespace std;

class Point;
class Box;
class Polytope;
class InterpolatedCurve;
class Zonotope;
class TaylorSet;
class GridCell;
class GridTreeSubset;
template<class BS> class ListSet;
template<class BS> class HybridBasicSet;
class DiscreteState;

    
//! \brief Class for plotting sets as a list of points.
class TextPlot 
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
    void draw(const GridCell&);
    void draw(const GridTreeSubset&);
    template<class BS> void draw(const ListSet<BS>&);    
    template<class BS> void draw(const HybridBasicSet<BS>&);
    template<class BS> void draw(const std::map< DiscreteState,BS >&);    
    void close();
  private:
    std::ofstream _fstream;
};

template<class BS> 
void TextPlot::draw(const ListSet<BS>& ls) {
    for(uint i=0; i!=ls.size(); ++i) {
        this->draw(ls[i]);
    }
}

template<class BS> 
void TextPlot::draw(const HybridBasicSet<BS>& hbs) {
    this->draw(hbs.continuous_state_set());
}

template<class DS> inline
void TextPlot::draw(const std::map<DiscreteState,DS>& hds) {
    for(typename std::map<DiscreteState,DS>::const_iterator loc_iter=hds.begin();
        loc_iter!=hds.end(); ++loc_iter) {
        this->draw(loc_iter->second);
    }
}

template<class SET> TextPlot& operator<<(TextPlot& g, const SET& set) { g.draw(set); return g; }

template<class SET> void textplot(const char* filename, const SET& set) { 
    TextPlot g(filename); g.draw(set); g.close(); }

} // namespace Ariadne

#endif // ARIADNE_TEXTPLOT_H
