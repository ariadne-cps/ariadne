/****************************************************************************
 *            epsstream.h
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_EPSSTREAM_H
#define ARIADNE_EPSSTREAM_H

/*! \file epsstream.h
 *  \brief Encapsulated postscript output.
 */
 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "base/stlio.h"
#include "linear_algebra/matrix.h"
#include "geometry/exceptions.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/box.h"
#include "geometry/rectangular_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedral_set.h"
#include "geometry/partition_tree_set.h"
#include "system/affine_map.h"
#include "output/colour.h"
#include "output/geometry2d.h"

namespace Ariadne {
  

    class LineStyle { bool _style; public: explicit LineStyle(bool ls) : _style(ls) { } operator bool() const { return this->_style; } };
    class FillColour : public Colour { public: FillColour(const Colour& fc) : Colour(fc) { } };

    inline LineStyle line_style(bool s) { return LineStyle(s); }
    inline FillColour fill_colour(const Colour& c) { return FillColour(c); }
    
    /*!\brief A stream for Encapsulated PostScript graphical output. */
    class epsstream
    {
     private:
      std::ostream* _os_ptr;
      
      static const uint xBBoxSide;
      static const uint yBBoxSide;
      static const double linewidth;
      static const double scale_dimension;
      
      PlanarProjectionMap p_map;
      Rectangle2d bbox;
     public:
      Colour line_colour;
      Colour fill_colour;
      bool line_style;
      bool fill_style;
     public:
      ~epsstream();
      epsstream();
      epsstream(std::ostream& os);
      void redirect(std::ostream& os);
      
      void write_header();
      void write_trailer();

      template<class R> void set_bounding_box(const Box<R>& bbox) {
        this->bbox=this->p_map(bbox); }

      void set_bounding_box(const Rectangle2d& bbox) {
        this->bbox=bbox; }

      void set_projection_map(const PlanarProjectionMap& p_map) {
        this->p_map=p_map; }
   
      
      void trace_scale(const char* x_name, const char* y_name,
                       const int& x_step=5, const int& y_step=5);
      
      const Rectangle2d& bounding_box() const { 
        return this->bbox; }

      const PlanarProjectionMap& projection_map() const { 
        return this->p_map; }

      std::ostream& ostream() { 
        return *this->_os_ptr;
      }

      void set_line_style(bool ls) { 
        this->line_style=ls;  }
    
      void set_fill_style(bool fs) {
        this->fill_style=fs; }
    
      void set_pen_colour(const Colour& pc) {
        this->line_colour=pc; }
    
      void set_fill_colour(const Colour& fc) {
        this->fill_colour=fc; }

      // FIXME: This is a hack to preserve the python interface
      void set_pen_colour(const char* pc) {
        this->line_colour=Colour(pc,0,0,0); }
    
      // FIXME: This is a hack to preserve the python interface
      void set_fill_colour(const char* fc) {
        this->fill_colour=Colour(fc,0,0,0); }

      void trace(const Point2d& pt);
      void trace(const Rectangle2d& r);
      void trace(const Zonotope2d& z);
      void trace(const Polygon2d& p);

      void draw(const Point2d& pt);
      void draw(const Segment2d& seg);
      void draw(const Rectangle2d& r);
      void draw(const Zonotope2d& z);
      void draw(const Polygon2d& p);

      void draw(std::vector<Rectangle2d>& rls);

      void fill();
      void stroke();
     private:
      static void _instantiate();
     public:     
      epsstream(const epsstream&); // no copy constructor
    };
    
  
    /*!\brief A stream for Encapsulated PostScript graphical output. */
    class epsfstream
      : public epsstream
    {
     private:
      std::ofstream* _ofs_ptr;
     public:
      ~epsfstream();
      epsfstream();
      
      template<class R> void open(const char* fn, const Box<R>& bbox);
      template<class R> void open(const char* fn, const Box<R>& bbox, unsigned int ix,  unsigned int iy);
      template<class R> void open(const char* fn, const Box<R>& bbox, const PlanarProjectionMap& p_map);
      void open(const char* fn, const Rectangle2d& bbox, const PlanarProjectionMap& p_map);
      void close();
     private:
      epsfstream(const epsfstream&); // no copy constructor
    };
 


    inline 
    epsstream& operator<<(epsstream& eps, std::ostream&(*f)(std::ostream&) ) {
      std::ostream& os=eps.ostream(); os << f; return eps; 
    }
      

    epsstream& operator<<(epsstream&, const LineStyle&); 
    epsstream& operator<<(epsstream&, const FillColour&);

    epsstream& operator<<(epsstream&, const char* s);

    template<class R> epsstream& operator<<(epsstream&, const Point<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const Segment<R>&);
    template<class R> epsstream& operator<<(epsstream&, const Box<R>&);
    template<class R> epsstream& operator<<(epsstream&, const Rectangle<R>&);
    template<class R> epsstream& operator<<(epsstream&, const Zonotope<R>&);
    template<class R> epsstream& operator<<(epsstream&, const Polytope<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const Polyhedron<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const RectangularSet<R>&);
    template<class R> epsstream& operator<<(epsstream&, const PolyhedralSet<R>&);
    template<class R> epsstream& operator<<(epsstream&, const BoxListSet<R>&); 
    template<class BS> epsstream& operator<<(epsstream&, const ListSet<BS>&); 
    template<class R> epsstream& operator<<(epsstream&, const GridCell<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const GridCellListSet<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const GridMaskSet<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const PartitionTreeSet<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const SetInterface<R>&); 

    template<class R> epsstream& operator<<(epsstream&, const FiniteGrid<R>&); 
    template<class R> epsstream& operator<<(epsstream&, const PartitionTree<R>&); 

    


    inline
    epsstream&
    operator<<(epsstream& eps, const LineStyle& ls)
    {
      eps.set_line_style(ls);
      return eps;
    }


    inline
    epsstream&
    operator<<(epsstream& eps, const FillColour& fc)
    {
      if(fc.transparant()) {
        eps.set_fill_style(false);
      } else {
        eps.set_fill_style(true);
        eps.set_fill_colour(fc);
      }
      return eps;
    }

  /*
    template<class R1, class R2>
    Point<R1> 
    approximate_point(const Point<R2>& pt) 
    {
      Point<R1> result(pt.dimension());
      for(size_type i=0; i!= result.dimension(); ++i) {
        result[i]=conv_approx<R1>(pt[i]);
      }
      return result;
    }
    
    template<class R1, class R2>
    PointList<R1> 
    approximate_point_list(const PointList<R2>& ptl) 
    {
      PointList<R1> result(ptl.dimension(),ptl.size());
      for(size_type i=0; i!= ptl.size(); ++i) {
        result.push_back(approximate_point<R1>(ptl[i]));
      }
      return result;
    }
  */
    
  

} // namespace Ariadne


#include "epsstream.inline.h"
#include "epsstream.template.h"

#endif /* ARIADNE_EPSSTREAM_H */
