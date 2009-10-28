/***************************************************************************
 *            graphics.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file graphics.h
 *  \brief Graphics class for drawting and outputting figures.
 */

#ifndef ARIADNE_GRAPHICS_H
#define ARIADNE_GRAPHICS_H

#include <iosfwd>
#include <string>
#include <vector>

#include "graphics_interface.h"

typedef unsigned int uint;

namespace Ariadne {

class Point;
class Box;
class Polytope;
class InterpolatedCurve;
class ProjectionFunction;

struct PlanarProjectionMap {
  public:
    PlanarProjectionMap(uint nn, uint ii, uint jj) : n(nn), i(ii), j(jj) { }
    uint argument_size() const { return this->n; }

    uint n, i, j;
};

inline std::ostream& operator<<(std::ostream& os, const PlanarProjectionMap& p) {
    return os << "PlanarProjectionMap( argument_size="<<p.argument_size()<<", x="<<p.i<<", y="<<p.j<<" )";
}


struct LineStyle { explicit LineStyle(bool ls) : _style(ls) { } operator bool() const { return this->_style; } private: bool _style; };
struct LineWidth { explicit LineWidth(double lw) : _width(lw) { } operator double() const { return this->_width; } private: double _width; };
struct LineColour : Colour { LineColour(const Colour& lc) : Colour(lc) { } LineColour(double r, double g, double b) : Colour(r,g,b) { } };
struct FillStyle { explicit FillStyle(bool fs) : _style(fs) { } operator bool() const { return this->_style; } private: bool _style; };
struct FillColour : Colour { FillColour(const Colour& fc) : Colour(fc) { } FillColour(double r, double g, double b) : Colour(r,g,b) { } };

inline LineStyle line_style(bool s) { return LineStyle(s); }
inline LineWidth line_width(double w) { return LineWidth(w); }
inline LineColour line_colour(const Colour& c) { return LineColour(c); }
inline LineColour line_colour(double r, double g, double b) { return LineColour(Colour(r,g,b)); }
inline FillStyle fill_style(bool s) { return FillStyle(s); }
inline FillColour fill_colour(const Colour& c) { return FillColour(c); }
inline FillColour fill_colour(double r, double g, double b) { return FillColour(Colour(r,g,b)); }
    

//! \brief Class for plotting figures.
class Figure 
    : public FigureInterface
{
  public:
    ~Figure();
    Figure();
    void set_projection_map(const ProjectionFunction&);
    void set_projection_map(const PlanarProjectionMap&);
    void set_bounding_box(const Box&);

    PlanarProjectionMap get_projection_map() const;
    Box get_bounding_box() const;

    void set_projection(uint, uint, uint);

    void set_x_axis_label(const string&);
    void set_y_axis_label(const string&);

    void set_line_style(bool);
    void set_line_width(double);
    void set_line_colour(Colour);
    void set_fill_style(bool);
    void set_fill_colour(Colour);

    void set_line_colour(double, double, double);
    void set_fill_colour(double, double, double);

    string get_x_axis_label() const;
    string get_y_axis_label() const;

    bool get_line_style() const;
    double get_line_width() const;
    Colour get_line_colour() const;
    bool get_fill_style() const;
    Colour get_fill_colour() const;

    void draw(const DrawableInterface& shape);

    void clear();
    void display();
    void write(const char* filename);
  public:
    class Data;
  public:
    void _paint_all(CanvasInterface& canvas); // Writes all shapes to the canvas
  private:
    Data* _data;
};


inline Figure& operator<<(Figure& g, const LineStyle& ls) { g.set_line_style(ls); return g; }
inline Figure& operator<<(Figure& g, const LineWidth& lw) { g.set_line_width(lw); return g; }
inline Figure& operator<<(Figure& g, const LineColour& lc) { g.set_line_colour(lc); return g; }
inline Figure& operator<<(Figure& g, const FillStyle& fs) { g.set_fill_style(fs); return g; }
inline Figure& operator<<(Figure& g, const FillColour& fc) { g.set_fill_colour(fc); return g; }

inline Figure& operator<<(Figure& fig, const DrawableInterface& shape) { fig.draw(shape); return fig; }

template<class SET> void plot(const char* filename, const SET& set) { 
    Figure g; draw(g,set); g.write(filename); }

template<class SET> void plot(const char* filename, const Colour& fc, const SET& set) { 
    Figure g; g.set_fill_colour(fc); draw(g,set); g.write(filename); }

template<class SET> void plot(const char* filename, const Box& bbox, const Colour& fc, const SET& set) { 
    Figure g; g.set_bounding_box(bbox); g.set_fill_colour(fc); draw(g,set); g.write(filename); }

template<class SET> void plot(const char* filename, const PlanarProjectionMap& pr, const Box& bbox, const Colour& fc, const SET& set) { 
    Figure g; g.set_projection_map(pr), g.set_bounding_box(bbox); g.set_fill_colour(fc); draw(g,set); g.write(filename); }

template<class SET1, class SET2> 
void plot(const char* filename, const Box& bbox, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2) { 
    Figure g; g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2); g.write(filename); }

template<class SET1, class SET2> 
void plot(const char* filename, const PlanarProjectionMap& pr, const Box& bbox, const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2) { 
    Figure g; g.set_bounding_box(bbox); g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2); g.write(filename); }

template<class SET1, class SET2, class SET3> 
void plot(const char* filename, const Box& bbox, 
          const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2, const Colour& fc3, const SET3& set3) 
{ 
    Figure g; g.set_bounding_box(bbox); 
    g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2); g.set_fill_colour(fc3); draw(g,set3); g.write(filename); 
}

template<class SET1, class SET2, class SET3> 
void plot(const char* filename, const PlanarProjectionMap& pr, const Box& bbox, 
          const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2, const Colour& fc3, const SET3& set3) 
{ 
    Figure g; g.set_projection_map(pr); g.set_bounding_box(bbox); 
    g.set_fill_colour(fc1); draw(g,set1); g.set_fill_colour(fc2); draw(g,set2); g.set_fill_colour(fc3); draw(g,set3); g.write(filename); 
}

template<class SET1, class SET2, class SET3, class SET4> 
void plot(const char* filename, const Box& bbox, 
          const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2, 
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4) 
{ 
    Figure g; g.set_bounding_box(bbox); 
    g.set_fill_colour(fc1); draw(g,set1); 
    g.set_fill_colour(fc2); draw(g,set2); 
    g.set_fill_colour(fc3); draw(g,set3); 
    g.set_fill_colour(fc4); draw(g,set4);
    g.write(filename); 
}

template<class SET1, class SET2, class SET3, class SET4> 
void plot(const char* filename, const PlanarProjectionMap& pr, const Box& bbox, 
          const Colour& fc1, const SET1& set1, const Colour& fc2, const SET2& set2, 
          const Colour& fc3, const SET3& set3, const Colour& fc4, const SET4& set4) 
{ 
    Figure g; g.set_projection_map(pr); g.set_bounding_box(bbox); 
    g.set_fill_colour(fc1); draw(g,set1); 
    g.set_fill_colour(fc2); draw(g,set2); 
    g.set_fill_colour(fc3); draw(g,set3); 
    g.set_fill_colour(fc4); draw(g,set4);
    g.write(filename); 
}

//  Luca's enhanced plot function.
//  if MAX_GRID_DEPTH is different from -1, plot the corresponding grid
/*
template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, const int& numVariables, const Box& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) { 
    // Assigns local variables
    Figure fig; 
    array<uint> xy(2,xaxis,yaxis);

    fig.set_projection_map(ProjectionFunction(xy,numVariables)); 
    fig.set_bounding_box(bbox); 
    
    fig.set_x_axis_label("x label");
    fig.set_y_axis_label("y label");
    
    // If the grid must be shown
    if (MAX_GRID_DEPTH >= 0)
    {
	// The rectangle to be drawn
	Box rect = Box(numVariables);
	// Chooses the fill colour
        fig << fill_colour(Colour(1.0,1.0,1.0));

	// Gets the number of times each variable interval would be divided by 2
        int numDivisions = MAX_GRID_DEPTH / numVariables;
	// Gets the step in the x direction, by 1/2^(numDivisions+h), where h is 1 if the step is to be further divided by 2, 0 otherwise
	double step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
	// Initiates the x position to the bounding box left bound
        double pos_x = bbox[0].lower();
        // Sets the rectangle 2-nd interval to the corresponding bounding box interval (while the >2 intervals are kept at [0,0])
	rect[yaxis] = bbox[1];
        // While between the interval
        while (pos_x < bbox[0].upper())
        {
	    rect[xaxis] = Interval(pos_x,pos_x+step_x); // Sets the rectangle x coordinate
	    pos_x += step_x; // Shifts the x position
	    fig << rect; // Appends the rectangle
        }

	// Repeats for the rectangles in the y direction
	double step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));  
        double pos_y = bbox[1].lower();
	rect[xaxis] = bbox[0];
        while (pos_y < bbox[1].upper())
        {
	    rect[yaxis] = Interval(pos_y,pos_y+step_y);
   	    fig << rect;
	    pos_y += step_y;
        }
    }
    // Draws and creates file
    fig.set_fill_colour(fc); 
    fig << set; 
    fig.write(filename); 
}
*/


} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
