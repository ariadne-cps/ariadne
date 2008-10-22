/****************************************************************************
 *            geometry2d.h
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

#ifndef ARIADNE_GEOMETRY2D_H
#define ARIADNE_GEOMETRY2D_H

/*! \file geometry2d.h
 *  \brief Classes for two-dimensional geometry output.
 */
 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "macros.h"
#include "point.h"

namespace Ariadne {

template<class X> class Vector;
class Point;
class Box;
class Zonotope;
class Polytope;
class InterpolatedCurve;
class ApproximateTaylorModel;

Polytope polytope(const Box& bx);
Polytope polytope(const Zonotope& z);
Polytope polytope(const Polytope& p);
Polytope polytope(const ApproximateTaylorModel& ts);

Point baricentre(const Polytope& p);
bool operator<(const Point& pt1, const Point& pt2);

Polytope& reduce2d(Polytope& p);
Float slope2d(const Point& pt1, const Point& pt2);

class Polytope
  : public std::vector<Point>  
{ 
 public: 
  Polytope() { }
  Polytope(uint d) { }
  Polytope(const std::vector<Point>& v) : std::vector<Point>(v) { }
  void new_vertex(const Point& v) { ARIADNE_ASSERT(this->size()==0 || v.dimension()==this->front().dimension()); this->push_back(v); }
  uint dimension() const { if(this->size()==0) { return 0; } return this->front().dimension(); }
  size_t number_of_vertices() const { return this->size(); }
  const std::vector<Point>& vertices() const { return *this; }
  std::vector<Point>& vertices() { return *this; }
  std::vector<Point>::const_iterator vertices_begin() const { return this->begin(); }
  std::vector<Point>::const_iterator vertices_end() const { return this->end(); }
  Box bounding_box() const;
};


class PlanarProjectionMap
{
 public:
  PlanarProjectionMap();
  PlanarProjectionMap(uint d, uint i, uint j);
  size_t argument_size() const;
  template<class X> Vector<X> operator() (const Vector<X>& v) const;
  Point operator() (const Point& pt) const;
  Box operator() (const Box& bx) const;
  Zonotope operator() (const Zonotope& z) const;
  Polytope operator() (const Polytope& z) const;
  InterpolatedCurve operator() (const InterpolatedCurve& c) const;
 private:
  friend std::ostream& operator<<(std::ostream&, const PlanarProjectionMap&);
 private:
  uint _d;
  uint _i;
  uint _j;
};
       


template<class X>
Vector<X>
PlanarProjectionMap::operator()(const Vector<X>& v) const 
{
  Vector<X> result(2); 
  ARIADNE_ASSERT(v.size()==this->_d);
  result[0]=v[this->_i]; 
  result[1]=v[this->_j]; 
  return result;
}



} // namespace Ariadne


#endif /* ARIADNE_OUTPUT_GEOMETRY2D_H */
