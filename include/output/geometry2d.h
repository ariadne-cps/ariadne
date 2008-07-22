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

#include "base/stlio.h"
#include "linear_algebra/matrix.h"
#include "geometry/exceptions.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/interpolated_curve.h"
#include "geometry/rectangular_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedral_set.h"
#include "geometry/partition_tree_set.h"
#include "system/affine_map.h"
#include "output/colour.h"


namespace Ariadne {
  

    struct Vector2d : array<double,2u> { Vector2d() : array<double,2u>(2u) { } };
    struct Point2d : array<double,2u> { Point2d() : array<double,2u>(2u) { } };
      

    std::ostream& operator<<(std::ostream&, const Point2d&);
    std::ostream& operator<<(std::ostream&, const Vector2d&);

    inline bool operator==(const Point2d& pt1, const Point2d& pt2) {
      return (pt1[0]==pt2[0]) &&  (pt1[1]==pt2[1]) ;
    }

    inline bool operator<(const Point2d& pt1, const Point2d& pt2) {
      return (pt1[0]<pt2[0]) || ( (pt1[0]==pt2[0]) && (pt1[1]<pt2[1]) );
    }

    inline Point2d& operator+=(Point2d& pt, const Vector2d& v) {
      pt[0]+=v[0]; pt[1]+=v[1]; return pt;
    }

    inline Point2d& operator-=(Point2d& pt, const Vector2d& v) {
      pt[0]-=v[0]; pt[1]-=v[1]; return pt;
    }

    class InterpolatedCurve2d {
     private:
      std::vector< Point2d > _points;
     public:
      typedef std::vector<Point2d>::const_iterator const_iterator;

      InterpolatedCurve2d(const Point2d& pt) : _points() { push_back(pt); }
      InterpolatedCurve2d(const Point2d& pt0, const Point2d& pt1) : _points() { push_back(pt0); push_back(pt1); }
      void push_back(const Point2d& pt) { _points.push_back(pt); }
      size_type size() const { return _points.size()-1; }
      const Point2d& operator[](dimension_type i) const { return _points[i]; }
      const_iterator begin() const { return _points.begin(); }
      const_iterator end() const { return --_points.end(); }
    };

    class Box2d {
     private:
      Point2d _lower_corner;
      Point2d _upper_corner;
     public:
      Box2d() { }
      Box2d(const Point2d& l, const Point2d& u) : _lower_corner(l), _upper_corner(u) { }
      dimension_type dimension() const { return 2; }
      const double& lower_bound(dimension_type i) const { return _lower_corner[i]; }
      const double& upper_bound(dimension_type i) const { return _upper_corner[i]; }
      double& lower_bound(dimension_type i) { return _lower_corner[i]; }
      double& upper_bound(dimension_type i) { return _upper_corner[i]; }
      void set_lower_bound(dimension_type i,const double& x) { _lower_corner[i]=x; }
      void set_upper_bound(dimension_type i,const double& x) { _upper_corner[i]=x; }
      friend bool operator==(const Box2d& r1, const Box2d& r2);
      friend bool operator<(const Box2d& r1, const Box2d& r2);
    };


    inline 
    bool operator==(const Box2d& r1, const Box2d& r2) {
      return r1._lower_corner==r2._lower_corner && r1._upper_corner==r2._upper_corner;
    }

    inline 
    bool operator<(const Box2d& r1, const Box2d& r2) {
      //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
      return r1._lower_corner<r2._lower_corner 
        || r1._lower_corner==r2._lower_corner && r1._upper_corner<r2._upper_corner; 
    }

    std::ostream& operator<<(std::ostream&, const Box2d&);



    struct Polygon2d {
      std::vector<Point2d> _vertices;

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


    struct Zonotope2d {
      Point2d centre;
      std::vector<Vector2d> generators;

      operator Polygon2d() const;
    };


    class PlanarProjectionMap;
    std::ostream& operator<<(std::ostream&, const PlanarProjectionMap&); 

    class PlanarProjectionMap
    {
     public:
      PlanarProjectionMap();
      PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j);
      template<class R> Vector2d operator() (const Vector<R>& v) const;
      template<class R> Point2d operator() (const Point<R>& pt) const;
      template<class R> InterpolatedCurve2d operator() (const InterpolatedCurve<R>& cv) const;
      template<class E> Box2d operator() (const BoxExpression<E>& re) const;
      template<class R> Zonotope2d operator() (const Zonotope<R>& z) const;
      template<class R> Polygon2d operator() (const Polytope<R>& p) const;
     private:
      friend std::ostream& operator<<(std::ostream&,const PlanarProjectionMap&);
     private:
      dimension_type _d;
      dimension_type _i;
      dimension_type _j;
    };
       


template<class R>
Vector2d
PlanarProjectionMap::operator()(const Vector<R>& v) const 
{
  Vector2d result; 
  ARIADNE_CHECK_SIZE(v,_d,"Point2d PlanarProjectionMap::operator()(Vector<R> v)");
  result[0]=approx<double>(v[_i]); 
  result[1]=approx<double>(v[_j]); 
  return result;
}

template<class R>
Point2d
PlanarProjectionMap::operator()(const Point<R>& pt) const 
{
  Point2d result; 
  ARIADNE_CHECK_DIMENSION(pt,_d,"Point2d PlanarProjectionMap::operator()(Point<R> pt)");
  result[0]=approx<double>(pt[_i]); 
  result[1]=approx<double>(pt[_j]); 
  return result;
}




template<class R> 
InterpolatedCurve2d
PlanarProjectionMap::operator()(const InterpolatedCurve<R>& curve) const
{
  const PlanarProjectionMap& self=*this;
  typename InterpolatedCurve<R>::const_iterator iter=curve.begin();
  typename InterpolatedCurve<R>::const_iterator end=curve.end();
  InterpolatedCurve2d result(self(iter->second));
  while(iter!=curve.end()) {
    ++iter;
    result.push_back(self(iter->second));
  }
  return result;
}

template<class E> 
Box2d 
PlanarProjectionMap::operator()(const BoxExpression<E>& re) const 
{
  Box2d result; 
  const E& r=re();
  ARIADNE_CHECK_DIMENSION(r,this->_d,"Box2d PlanarProjectionMap::operator()(Box<R> r)");
  result.lower_bound(0)=approx<double>(r.lower_bound(this->_i));
  result.upper_bound(0)=approx<double>(r.upper_bound(this->_i));
  result.lower_bound(1)=approx<double>(r.lower_bound(this->_j));
  result.upper_bound(1)=approx<double>(r.upper_bound(this->_j));
  return result;
}


template<class R> 
Zonotope2d 
PlanarProjectionMap::operator() (const Zonotope<R>& z) const 
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  const PlanarProjectionMap& map=*this;
  Zonotope2d result;
  result.centre=map(z.centre());
  for(size_type i=0; i!=z.number_of_generators(); ++i) {
    Vector<R> v=z.generators().column(i);
    result.generators.push_back(map(v));
  }
  return result;  
}


template<class R> 
Polygon2d 
PlanarProjectionMap::operator() (const Polytope<R>& p) const 
{
  Polygon2d result;
  for(size_type i=0; i!=p.number_of_vertices(); ++i) {
    Point<R> v=p.vertex(i);
    Point2d pt=(*this)(v);
    result.new_vertex(pt);
  }
  result.reduce();
  return result;
}


} // namespace Ariadne


#endif /* ARIADNE_OUTPUT_GEOMETRY2D_H */
