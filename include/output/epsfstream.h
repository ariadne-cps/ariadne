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

#ifndef ARIADNE_EPSFSTREAM_H
#define ARIADNE_EPSFSTREAM_H

/*! \file epsfstream.h
 *  \brief Encapsulated postscript output.
 */
 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "../base/stlio.h"
#include "../numeric/conversion.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/exceptions.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/rectangular_set.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"
#include "../geometry/polyhedral_set.h"
#include "../geometry/partition_tree_set.h"
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
    bool operator==(const Point2d& pt1, const Point2d& pt2) {
      return pt1[0]==pt2[0] && pt1[1]==pt2[1];
    }

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
      friend bool operator==(const Rectangle2d& r1, const Rectangle2d& r2);
      friend bool operator<(const Rectangle2d& r1, const Rectangle2d& r2);
    };


    inline 
    bool operator==(const Rectangle2d& r1, const Rectangle2d& r2) {
      return r1._lower_corner==r2._lower_corner && r1._upper_corner==r2._upper_corner;
    }

    inline 
    bool operator<(const Rectangle2d& r1, const Rectangle2d& r2) {
      return r1._lower_corner<r2._lower_corner 
        || r1._lower_corner==r2._lower_corner && r1._upper_corner<r2._upper_corner; 
    }

    std::ostream& operator<<(std::ostream&, const Rectangle2d&);


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


    template<class R1, class R2>
    Geometry::Point<R1> 
    approximate_point(const Geometry::Point<R2>& pt) 
    {
      Geometry::Point<R1> result(pt.dimension());
      for(size_type i=0; i!= result.dimension(); ++i) {
        result[i]=Numeric::conv_approx<R1>(pt[i]);
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
    

    class PlanarProjectionMap;
    std::ostream& operator<<(std::ostream&, const PlanarProjectionMap&); 

    class PlanarProjectionMap
    {
     public:
      PlanarProjectionMap() : _d(2), _i(0), _j(1) { }
      PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j)
        : _d(d), _i(i), _j(j) 
      { 
        if(i>=d || j>=d) { 
          using Geometry::InvalidCoordinate; 
          ARIADNE_THROW(InvalidCoordinate,"PlanarProjectionMap::PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j)",
                        "d="<<d<<", i="<<i<<", j="<<j); }
      }
      template<class R> Point2d operator() (const Geometry::Point<R>& pt) const {
        Point2d result(2); 
        ARIADNE_CHECK_DIMENSION(pt,_d,"Point2d PlanarProjectionMap::operator()(Point<R> pt)");
        result[0]=Numeric::conv_approx<double>(pt[_i]); 
        result[1]=Numeric::conv_approx<double>(pt[_j]); 
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
      template<class E> Rectangle2d operator() (const Geometry::RectangleExpression<E>& re) const {
        Rectangle2d result(2); 
        const E& r=re();
        ARIADNE_CHECK_DIMENSION(r,this->_d,"Rectangle2d PlanarProjectionMap::operator()(Rectangle<R> r)");
        result.lower_bound(0)=Numeric::conv_approx<double>(r.lower_bound(this->_i));
        result.upper_bound(0)=Numeric::conv_approx<double>(r.upper_bound(this->_i));
        result.lower_bound(1)=Numeric::conv_approx<double>(r.lower_bound(this->_j));
        result.upper_bound(1)=Numeric::conv_approx<double>(r.upper_bound(this->_j));
        return result;
      }
      template<class RC,class RG> Polygon2d operator() (const Geometry::Zonotope<RC,RG>& z) const {
        Geometry::Point<Numeric::Rational> c=approximate_value(z.centre());
        LinearAlgebra::Matrix<Numeric::Rational> G=approximate_value(z.generators());
        return this->operator()(Geometry::Zonotope<Numeric::Rational>(c,G).vertices());
      }
      template<class R> Polygon2d operator() (const Geometry::Polytope<R>& p) const {
        return this->operator()(p.vertices());
      }
     private:
      friend std::ostream& operator<<(std::ostream&,const PlanarProjectionMap&);
     private:
      dimension_type _d;
      dimension_type _i;
      dimension_type _j;
    };
         

    Point2d baricentre(const Polygon2d& vertices);

    class LineWidth {
     public:
      explicit LineWidth(double x) : _value(x) { }
      double value() const { return this->_value; }
     private:
      double _value;
    };

    LineWidth line_width(double x);
  
    class FillColour {
     public:
      explicit FillColour(std::string x) : _value(x) { }
      std::string value() const { return this->_value; }
     private:
      std::string _value;
    };

    FillColour fill_colour(double x);
        

    
    class epsfstream
      : private std::ofstream 
    {
     private:
      static const uint xBBoxSide;
      static const uint yBBoxSide;
      static const double linewidth;
      
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
      void open(const char* fn, const Geometry::Rectangle<R>& bbox);
      
      template<class R>
      void open(const char* fn, const Geometry::Rectangle<R>& bbox, 
                unsigned int ix,  unsigned int iy);
      
      template<class R>
      void open(const char* fn, const Geometry::Rectangle<R>& bbox, 
                const PlanarProjectionMap& p_map);
      
      void open(const char* fn, const Rectangle2d& bbox,
                const PlanarProjectionMap& p_map);

      void close();
      
      void trace_scale(const char* x_name, const char* y_name,
                       const int& x_step=5, const int& y_step=5);
      
      const Rectangle2d& bounding_box() const { 
        return this->bbox; }

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

      void trace(const Point2d& pt);
      void trace(const Rectangle2d& r);
      void trace(const Polygon2d& p);

      void draw(const Point2d& pt);
      void draw(const Rectangle2d& r);
      void draw(const Polygon2d& p);

      void draw(std::vector<Rectangle2d>& rls);

      void fill();
      void stroke();
      
      double(*fff)(void) ;
      double(*ffz) ;
     private:
      friend epsfstream& operator<<(epsfstream&, std::ostream&(*)(std::ostream&) );
     private:
      void write_header();
      void write_trailer();
      epsfstream(const epsfstream&); // no copy constructor
    };
    
    inline 
    epsfstream& operator<<(epsfstream& eps, std::ostream&(*f)(std::ostream&) ) {
      std::ostream& os(eps); os << f; return eps; 
    }
      

    epsfstream& operator<<(epsfstream&, const LineWidth&); 
    epsfstream& operator<<(epsfstream&, const FillColour&);

    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::Point<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::Rectangle<R>&);
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::Zonotope<R,R>&);
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::Zonotope<Numeric::Interval<R>,R>&);
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::Zonotope< Numeric::Interval<R>,Numeric::Interval<R> >&);
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::Polytope<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::Polyhedron<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::RectangularSet<R>&);
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::PolyhedralSet<R>&);
    template<class BS> epsfstream& operator<<(epsfstream&, const Geometry::ListSet<BS>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::GridCellListSet<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::GridMaskSet<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::PartitionTreeSet<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::SetInterface<R>&); 

    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::FiniteGrid<R>&); 
    template<class R> epsfstream& operator<<(epsfstream&, const Geometry::PartitionTree<R>&); 

    

    template<class R>
    void 
    epsfstream::open(const char* fn, const Geometry::Rectangle<R>& bbox)
    {
      PlanarProjectionMap p_map(bbox.dimension(),0,1);
      this->open(fn,p_map(bbox),p_map);
    }

    template<class R>
    void 
    epsfstream::open(const char* fn, const Geometry::Rectangle<R>& bbox,
                     unsigned int ix,  unsigned int iy)
    {
      PlanarProjectionMap p_map(bbox.dimension(),ix,iy);
      this->open(fn,p_map(bbox),p_map);
    }

    template<class R>
    void 
    epsfstream::open(const char* fn, const Geometry::Rectangle<R>& bbox,
                     const PlanarProjectionMap& p_map)
    {
      this->open(fn,p_map(bbox),p_map);
    }


    inline 
    FillColour 
    fill_colour(const std::string& x) 
    {
      return FillColour(x);
    }

    inline 
    LineWidth 
    line_width(double x) 
    {
      return LineWidth(x);
    }

    inline 
    LineWidth 
    line_style(bool b) 
    {
      return LineWidth(double(b));
    }

    inline
    epsfstream&
    operator<<(epsfstream& eps, const LineWidth& lw)
    {
      if(lw.value()==0.0) {
        eps.set_line_style(false);
      } else {
        eps.set_line_style(true);
        //eps.set_line_width(lw.value());
      }
      return eps;
    }


    inline
    epsfstream&
    operator<<(epsfstream& eps, const FillColour& fc)
    {
      if(fc.value()=="none" || fc.value()=="transparant") {
        eps.set_fill_style(false);
      } else {
        eps.set_fill_style(true);
        eps.set_fill_colour(fc.value().c_str());
      }
      return eps;
    }


    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Point<R>& pt) 
    {
      Point2d dpt=eps.projection_map()(pt);
      eps.draw(dpt);
      return eps;
    }

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Rectangle<R>& r) 
    {
      Rectangle2d dr=eps.projection_map()(r);
      eps.draw(dr);
      return eps;
    }
    
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::RectangularSet<R>& rs)
    {
      return eps << Geometry::Rectangle<R>(rs);
    }

    
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Zonotope< Numeric::Interval<R>,Numeric::Interval<R> >& iz)
    { 
      Geometry::Zonotope<Numeric::Interval<R>,R> ez=Geometry::over_approximation(iz);
      Geometry::Zonotope<Numeric::Rational> qz=Geometry::over_approximation(ez);
      Polygon2d vertices=eps.projection_map()(qz);      
      eps.draw(vertices);
      return eps;
   }
       
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Zonotope<Numeric::Interval<R>,R>& ez)
    { 
      Geometry::Zonotope<R> z=Geometry::over_approximation(ez);
      Geometry::Zonotope<Numeric::Rational> qz(z);
      Polygon2d vertices=eps.projection_map()(qz);      
      eps.draw(vertices);
      return eps;
    }
       
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Zonotope<R>& z)
    { 
      Polygon2d vertices=eps.projection_map()(Geometry::Zonotope<Numeric::Rational>(z));      
      eps.draw(vertices);
      return eps;
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
    operator<<(epsfstream& eps, const Geometry::Polytope<R>& p)
    {
      eps.draw(eps.projection_map()(p.vertices()));
      return eps;
    }
    
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Polyhedron<R>& p)
    {
      return eps << Geometry::Polytope<Numeric::Rational>(p);
    }

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::PolyhedralSet<R>& ps)
    {
      return eps << Geometry::Polyhedron<R>(ps);
    }

    
    template<class BS> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::ListSet<BS>& ds)
    {
      typedef typename Geometry::ListSet<BS>::const_iterator const_iterator;
      if(eps.fill_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          eps.trace(eps.projection_map()(*set_iter));
          eps.fill();
        }
      }
      if(eps.line_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          eps.trace(eps.projection_map()(*set_iter));
          eps.stroke();
        }
      }
      return eps;
    }
    

 
    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::GridCell<R>& bs)
    {
      return eps << Geometry::Rectangle<R>(bs);
    }
    

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::GridBlock<R>& bs)
    {
      return eps << Geometry::Rectangle<R>(bs);
    }
    

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::GridCellListSet<R>& ds)
    {
      std::vector<Rectangle2d> rl(ds.size());
      std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
      eps.draw(rl);
      return eps;
    }
    

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::GridMaskSet<R>& ds)
    {
      std::vector<Rectangle2d> rl(ds.size());
      std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
      eps.draw(rl);
      return eps;
    }
    
 

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::PartitionTreeSet<R>& ds)
    {
      std::vector<Rectangle2d> rl(ds.size());
      std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
      eps.draw(rl);
      return eps;
    }

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::SetInterface<R>& set)
    {
      using namespace Geometry;
      typedef Numeric::Interval<R> I;
      
      if(dynamic_cast<const RectangularSet<R>*>(&set)) {
        return eps << dynamic_cast<const RectangularSet<R>&>(set);
      } else if(dynamic_cast<const PolyhedralSet<R>*>(&set)) {
        return eps << dynamic_cast<const PolyhedralSet<R>&>(set);
      } else if(dynamic_cast<const ListSet< Rectangle<R> >*>(&set)) {
        return eps << dynamic_cast<const ListSet< Rectangle<R> >&>(set);
      } else if(dynamic_cast<const ListSet< Zonotope<R,R> >*>(&set)) {
        return eps << dynamic_cast<const ListSet< Zonotope<R,R> >&>(set);
      } else if(dynamic_cast<const ListSet< Zonotope<I,R> >*>(&set)) {
        return eps << dynamic_cast<const ListSet< Zonotope<I,R> >&>(set);
      } else if(dynamic_cast<const ListSet< Zonotope<I,I> >*>(&set)) {
        return eps << dynamic_cast<const ListSet< Zonotope<I,I> >&>(set);
      } else if(dynamic_cast<const GridCellListSet<R>*>(&set)) {
        return eps << dynamic_cast<const GridCellListSet<R>&>(set);
      } else if(dynamic_cast<const GridMaskSet<R>*>(&set)) {
        return eps << dynamic_cast<const GridMaskSet<R>&>(set);
      } else if(dynamic_cast<const PartitionTreeSet<R>*>(&set)) {
        return eps << dynamic_cast<const PartitionTreeSet<R>&>(set);
      }  else {
        Rectangle<R> bb;
        try {
          bb=set.bounding_box();
        } 
        catch(Geometry::UnboundedSet& e) {
          if(set.dimension()==2) {
            Rectangle2d bbox=eps.bounding_box();
            bb=Geometry::Rectangle<R>(2);
            bb.set_lower_bound(0,bbox.lower_bound(0));
            bb.set_upper_bound(0,bbox.upper_bound(0));
            bb.set_lower_bound(1,bbox.lower_bound(1));
            bb.set_upper_bound(1,bbox.upper_bound(1));
          } else {
            throw e;
          }
        }
        Geometry::PartitionScheme<R> ps(bb);
        int depth=16;
        Geometry::PartitionTreeSet<R> pts=Geometry::outer_approximation(set,ps,depth);
        return eps << pts;
      }
    }

    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::FiniteGrid<R>& fg)
    {
      bool fill_style=eps.fill_style;
      if(fill_style) { eps.fill_style=false; }
      Geometry::GridCellListSet<R> gcls(fg.grid());
      gcls.adjoin(Geometry::GridBlock<R>(fg.grid(),fg.lattice_block()));
      eps << gcls;
      if(fill_style) { eps.fill_style=true; }
      return eps;
    }


    template<class R> inline
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::PartitionTree<R>& pt)
    {
      bool fill_style=eps.fill_style;
      if(fill_style) { eps.fill_style=false; }
      for(typename Geometry::PartitionTree<R>::const_iterator iter = pt.begin(); iter!=pt.end(); ++iter) {
        eps << Geometry::Rectangle<R>(*iter);
      }
      if(fill_style) { eps.fill_style=true; }
      return eps;
    }
    
  }
}


#endif /* ARIADNE_EPSFSTREAM_H */
