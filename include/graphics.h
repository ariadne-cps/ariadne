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
 *  \brief Graphics class for drawing and outputting figures.
 */

#ifndef ARIADNE_GRAPHICS_H
#define ARIADNE_GRAPHICS_H

#include <cstdarg>
#include <cstdio>
#include <iosfwd>
#include <string>
#include <vector>

#include "graphics_interface.h"
#include "discrete_state.h"
#include "box.h"

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

//! \brief Plots figures for each couple of coordinates, with graphic box automatically chosen as the bounding box of the set. Also, saves the figures on a given folder under the current path.
template<class SET> 
void plot(const string& foldername, const string& filename, const SET& set)  
{
	// If a set exists
	if (set.size()>0)
	{
		// Gets the bounds of the set
		std::map<DiscreteState,Box> set_bounds = set.bounding_box();

		// Gets the number of variables (NOTE: assumed equal for each mode)
		uint numvar = set_bounds.begin()->second.dimension();	

		// For each variable (last excluded)
		for (uint x=0;x<numvar-1;x++)
		{
			// For each following variable
			for (uint y=x+1;y<numvar;y++)
			{
				// Initializes the iterator
				std::map<DiscreteState,Box>::const_iterator it = set_bounds.begin();

				// Sets the initial value for the graphics box	
				Box graphics_box(2);
				graphics_box[0] = it->second[x];
				graphics_box[1] = it->second[y];

				// For each other location			
				while (++it != set_bounds.end())
				{
					// Gets the bounding box for the location
					Box evaluation_box(2);
					evaluation_box[0] = it->second[x]; 
					evaluation_box[1] = it->second[y];
			
					// Enlarges the graphics box

					// For each axis
					for (uint u=0;u<2;u++)
					{
						// If the evaluated box has higher upper bound, extends the graphics box
						if (evaluation_box[u].upper() > graphics_box[u].upper())
							graphics_box[u].set_upper(evaluation_box[u].upper());
						// If the evaluated box has lower lower bound, extends the graphics box
						if (evaluation_box[u].lower() < graphics_box[u].lower())
							graphics_box[u].set_lower(evaluation_box[u].lower());
					}
				}
	
				// Plots the result

				// Assigns local variables
				Figure fig; 
				array<uint> xy(2,x,y);

				fig.set_projection_map(ProjectionFunction(xy,numvar)); 
				fig.set_bounding_box(graphics_box); 

				// Appends the set, with the desired fill color
				fig.set_fill_colour(Colour(1.0,0.75,0.0)); 
				draw(fig,set); 

				// If there are more than two variables, prints the variable numbers				
				char num_char[6] = "";
				if (numvar>2)
					sprintf(num_char,"[%u,%u]",x,y);
				// Writes the figure file
				fig.write((foldername+"/"+filename+num_char).c_str()); 
			}
		} 
	}
	else
		ARIADNE_WARN("Empty set, no plotting produced.");
}


} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
