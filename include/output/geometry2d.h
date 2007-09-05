/****************************************************************************
 *            output/geometry2d.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_OUTPUT_GEOMETRY2D_H
#define ARIADNE_OUTPUT_GEOMETRY2D_H

/*! \file output/geometry2d.h
 *  \brief Classes for two-dimensional geometry output.
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
#include "../output/colour.h"

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
      Point2d baricentre() const;
      Polygon2d& reduce();
     private:
      double slope(const Point2d& pt1, const Point2d& pt2) {
        return (pt2[1]-pt1[1])/(pt2[0]-pt1[0]); }
    };



    class PlanarProjectionMap;
    std::ostream& operator<<(std::ostream&, const PlanarProjectionMap&); 

    class PlanarProjectionMap
    {
     public:
      PlanarProjectionMap();
      PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j);
      template<class R> Point2d operator() (const Geometry::Point<R>& pt) const;
      template<class R> Polygon2d operator() (const Geometry::PointList<R>& ptl) const;
      template<class E> Rectangle2d operator() (const Geometry::RectangleExpression<E>& re) const;
      template<class RC,class RG> Polygon2d operator() (const Geometry::Zonotope<RC,RG>& z) const;
      template<class R> Polygon2d operator() (const Geometry::Polytope<R>& p) const;
     private:
      friend std::ostream& operator<<(std::ostream&,const PlanarProjectionMap&);
     private:
      dimension_type _d;
      dimension_type _i;
      dimension_type _j;
    };
       
  }
}


namespace Ariadne {

template<class R>
Output::Point2d
Output::PlanarProjectionMap::operator()(const Geometry::Point<R>& pt) const 
{
  Point2d result(2); 
  ARIADNE_CHECK_DIMENSION(pt,_d,"Point2d PlanarProjectionMap::operator()(Point<R> pt)");
  result[0]=Numeric::conv_approx<double>(pt[_i]); 
  result[1]=Numeric::conv_approx<double>(pt[_j]); 
  return result;
}


template<class R> 
Output::Polygon2d 
Output::PlanarProjectionMap::operator()(const Geometry::PointList<R>& ptl) const
{
  Polygon2d result(2);
  for(size_type i=0; i!=ptl.size(); ++i) { 
    result.new_vertex(this->operator()(ptl[i]));
  }
  result.reduce();
  return result;
}


template<class E> 
Output::Rectangle2d 
Output::PlanarProjectionMap::operator()(const Geometry::RectangleExpression<E>& re) const 
{
  Rectangle2d result(2); 
  const E& r=re();
  ARIADNE_CHECK_DIMENSION(r,this->_d,"Rectangle2d PlanarProjectionMap::operator()(Rectangle<R> r)");
  result.lower_bound(0)=Numeric::conv_approx<double>(r.lower_bound(this->_i));
  result.upper_bound(0)=Numeric::conv_approx<double>(r.upper_bound(this->_i));
  result.lower_bound(1)=Numeric::conv_approx<double>(r.lower_bound(this->_j));
  result.upper_bound(1)=Numeric::conv_approx<double>(r.upper_bound(this->_j));
  return result;
}


template<class RC,class RG> 
Output::Polygon2d 
Output::PlanarProjectionMap::operator() (const Geometry::Zonotope<RC,RG>& z) const 
{
  Geometry::Point<Numeric::Rational> c=Geometry::approximation(z.centre());
  LinearAlgebra::Matrix<Numeric::Rational> G=LinearAlgebra::approximation<Numeric::Rational>(z.generators());
  return this->operator()(Geometry::Zonotope<Numeric::Rational>(c,G).vertices());
}


template<class R> 
Output::Polygon2d 
Output::PlanarProjectionMap::operator() (const Geometry::Polytope<R>& p) const 
{
  return this->operator()(p.vertices());
}


}


#endif /* ARIADNE_OUTPUT_GEOMETRY2D_H */
