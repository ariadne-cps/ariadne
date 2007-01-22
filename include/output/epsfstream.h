/****************************************************************************
 *            epsfstream.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

/*! \file epsfstream.h
 *  \brief Encapsulated postscript output.
 */
 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "../base/stlio.h"
#include "../numeric/conversion.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"
#include "../geometry/partition_tree_set.h"
#include "../linear_algebra/matrix.h"
#include "../system/affine_map.h"

#define SCALE_DIMENSION 3.5 

// FIXME: This should not be necessary
namespace Ariadne {
  namespace Numeric {

    template<> inline double min_exact(const double& x1, const double& x2) {
      return std::min(x1,x2);
    }
    
    template<> inline double max_exact(const double& x1, const double& x2) {
      return std::max(x1,x2);
    }
    
  }
}

namespace Ariadne {
  namespace Output {
    
    class Point2d {
     private:
      double _coordinates[2];
     public:
      Point2d(dimension_type d=2) { assert(d==2); }
      Point2d(double x, double y) { _coordinates[0]=x; _coordinates[1]=y; }
      dimension_type dimension() const { return 2; }
      const double& operator[](dimension_type i) const { return _coordinates[i]; }
      double& operator[](dimension_type i) { return _coordinates[i]; }
    };

    std::ostream& operator<<(std::ostream&, const Point2d&);

    inline 
    bool operator<(const Point2d& pt1, const Point2d& pt2) {
      return (pt1[0]<pt2[0]) || ( (pt1[0]==pt2[0]) && (pt1[1]<pt2[1]) );
    }

    class Rectangle2d {
     private:
      Point2d _lower_corner;
      Point2d _upper_corner;
     public:
      Rectangle2d(dimension_type d=2) { assert(d==2); }
      Rectangle2d(const Point2d& l, const Point2d& u) : _lower_corner(l), _upper_corner(u) { }
      dimension_type dimension() const { return 2; }
      const double& lower_bound(dimension_type i) const { return _lower_corner[i]; }
      const double& upper_bound(dimension_type i) const { return _upper_corner[i]; }
      double& lower_bound(dimension_type i) { return _lower_corner[i]; }
      double& upper_bound(dimension_type i) { return _upper_corner[i]; }
      void set_lower_bound(dimension_type i,const double& x) { _lower_corner[i]=x; }
      void set_upper_bound(dimension_type i,const double& x) { _upper_corner[i]=x; }
    };

    class Polygon2d {
     private:
      std::vector<Point2d> _vertices;
     public:
      Polygon2d(dimension_type d=2) { assert(d==2); }
      Polygon2d(dimension_type d,size_type nv) { assert(d==2); _vertices.resize(nv); }
      dimension_type dimension() const { return 2; }
      size_type size() const { return _vertices.size(); }
      const Point2d& operator[](dimension_type i) const { return _vertices[i]; }
      Point2d& operator[](dimension_type i) { return _vertices[i]; }
      void new_vertex(const Point2d& pt) { _vertices.push_back(pt); }
      Polygon2d& reduce();
     private:
      double slope(const Point2d& pt1, const Point2d& pt2) {
        return (pt2[1]-pt1[1])/(pt2[0]-pt1[0]); }
    };


  typedef struct{
      size_t pos;
      double radiant;
    } radiant_pointer_type;


    inline bool
    is_smaller_than(const radiant_pointer_type &a, 
                    const radiant_pointer_type &b) {
      return (a.radiant<b.radiant);
    }
   

    template<class R1, class R2>
    Geometry::Point<R1> 
    approximate_point(const Geometry::Point<R2>& pt) 
    {
      Geometry::Point<R1> result(pt.dimension());
      for(size_type i=0; i!= result.dimension(); ++i) {
        result[i]=conv_approx<R1>(pt[i]);
      }
      return result;
    }
    
    template<class R1, class R2>
    Geometry::PointList<R1> 
    approximate_point_list(const Geometry::PointList<R2>& ptl) 
    {
      Geometry::PointList<R1> result(ptl.dimension(),ptl.size());
      for(size_type i=0; i!= ptl.size(); ++i) {
        result.push_back(approximate_point<R1>(ptl[i]));
      }
      return result;
    }
    


    class PlanarProjectionMap
    {
     public:
      PlanarProjectionMap() : _d(2), _i(0), _j(1) { }
      PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j)
        : _d(d), _i(i), _j(j) { if(i>=d || j>=d) { throw Geometry::InvalidCoordinate(__PRETTY_FUNCTION__); } }
      template<class R> Point2d operator() (const Geometry::Point<R>& pt) const {
        Point2d result(2); 
        check_dimension(pt,_d,__PRETTY_FUNCTION__);
        result[0]=conv_approx<double>(pt[_i]); 
        result[1]=conv_approx<double>(pt[_j]); 
        return result;
      }
      template<class R> Polygon2d operator() (const Geometry::PointList<R>& ptl) const {
        Polygon2d result(2);
        for(size_type i=0; i!=ptl.size(); ++i) { 
          result.new_vertex(this->operator()(ptl[i]));
        }
        result.reduce();
        return result;
      }
      template<class R> Rectangle2d operator() (const Geometry::Rectangle<R>& r) const {
        Rectangle2d result(2); 
        check_dimension(r,this->_d);
        result.lower_bound(0)=conv_approx<double>(r.lower_bound(0));
        result.upper_bound(0)=conv_approx<double>(r.upper_bound(0));
        result.lower_bound(1)=conv_approx<double>(r.lower_bound(1));
        result.upper_bound(1)=conv_approx<double>(r.upper_bound(1));
        return result;
      }
      template<class R> Polygon2d operator() (const Geometry::Zonotope<R>& z) const {
        return this->operator()(Geometry::Zonotope<Rational>(z).vertices());
      }
      template<class R> Polygon2d operator() (const Geometry::Polytope<R>& p) const {
        return this->operator()(p.vertices());
      }
     private:
      dimension_type _d;
      dimension_type _i;
      dimension_type _j;
      };
    
    Point2d baricentre(const Polygon2d& vertices);


    Polygon2d
    order_around_a_point(const Polygon2d& vertices, 
                         const Point2d &centre);

    
    
    class epsfstream
      : public std::ofstream 
    {
     private:
      static const uint xBBoxSide=300;
      static const uint yBBoxSide=300;
      static const double linewidth=0.0000001;
      
      PlanarProjectionMap p_map;
      Rectangle2d bbox;
     public:
      std::string line_colour;
      std::string fill_colour;
      bool line_style;
      bool fill_style;
     public:
      ~epsfstream();

      epsfstream();
      
      template<class R>
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox, 
                 const unsigned int &ix=0,  const unsigned int& iy=1);
      
      template<class R>
      void open(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox, 
                const unsigned int &ix=0,  const unsigned int& iy=1);
      
      void open(const char* fn, const Rectangle2d& bbox,
                const PlanarProjectionMap& p_map);

      void close();
      
      void trace_scale(const char* x_name, const char* y_name,
                       const int& x_step=5, const int& y_step=5);
      
      const PlanarProjectionMap& projection_map() const { 
        return this->p_map; }

      void set_line_style(bool ls) { 
        this->line_style=ls;  }
    
      void set_fill_style(bool fs) {
        this->fill_style=fs; }
    
      void set_pen_colour(const char* pc) {
        this->line_colour=pc; }
    
      void set_fill_colour(const char* fc) {
        this->fill_colour=fc; }
     private:
      void write_header();
      void write_trailer();
      epsfstream(const epsfstream&); // no copy constructor
    };
    
    epsfstream& trace(epsfstream& eps, const Point2d& pt);
    epsfstream& trace(epsfstream& eps, const Rectangle2d& r);
    epsfstream& trace(epsfstream& eps, const Polygon2d& vertices);

    epsfstream& draw(epsfstream& eps, const Point2d& pt);
    epsfstream& draw(epsfstream& eps, const Rectangle2d& r);
    epsfstream& draw(epsfstream& eps, const Polygon2d& vertices);
    
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::Point<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::Rectangle<R>&);
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::Zonotope<R>&);
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::Polytope<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::Polyhedron<R>&); 
    template<class R,template<class> class BS> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::ListSet<R,BS>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::GridCellListSet<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::GridMaskSet<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::PartitionTree<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Ariadne::Geometry::PartitionTreeSet<R>&); 

    

    template<class R>
    epsfstream::epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox, 
                           const unsigned int &ix,  const unsigned int& iy)
      : std::ofstream(), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
    {
      this->open(fn,bbox,ix,iy);
    }


    template<class R>
    void 
    epsfstream::open(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox,
                     const unsigned int &ix,  const unsigned int& iy)
    {
      PlanarProjectionMap p_map(bbox.dimension(),ix,iy);
      this->open(fn,p_map(bbox),p_map);
    }



    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Point<R>& pt) 
    {
      return draw(eps, eps.projection_map()(pt));
    }

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Rectangle<R>& r) 
    {
      Rectangle2d dr=eps.projection_map()(r);
      return draw(eps,dr);
    }
    
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Zonotope<R>& z)
    {
      Polygon2d vertices=eps.projection_map()(Geometry::Zonotope<Rational>(z).vertices());      
      return draw(eps,vertices);
    }
       
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Parallelotope<R>& p)
    {
      const Geometry::Zonotope<R>& z=p;
      return eps << z;
    }
       
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Polytope<R>& p)
    {
      return draw(eps,eps.projection_map()(p.vertices()));      
    }
    
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Polyhedron<R>& p)
    {
      return draw(eps,eps.projection_map()(Geometry::Polytope<Rational>(Geometry::Polyhedron<Rational>(p))));      
    }
    
    template<class R, template<class> class BS> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::ListSet<R,BS>& ds)
    {
      typedef typename Ariadne::Geometry::ListSet<R,BS>::const_iterator const_iterator;
      if(eps.fill_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          trace(eps,eps.projection_map()(*set_iter)) << eps.fill_colour << " fill\n";
        }
      }
      if(eps.line_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          trace(eps,eps.projection_map()(*set_iter)) << eps.line_colour << " stroke\n";
        }
      }
      return eps;
    }
    

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::GridMaskSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    
 
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::GridCellListSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::PartitionTreeSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::PartitionTree<R>& pt)
    {
      for(typename Ariadne::Geometry::PartitionTree<R>::const_iterator iter = pt.begin(); iter!=pt.end(); ++iter) {
        eps << Ariadne::Geometry::Rectangle<R>(*iter);
      }
      return eps;
    }
    
  }
}
