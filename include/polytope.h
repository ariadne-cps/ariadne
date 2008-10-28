/****************************************************************************
 *            polytope.h
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

#ifndef ARIADNE_POLYTOPE_H
#define ARIADNE_POLYTOPE_H

/*! \file polytope.h
 *  \brief Polytope class for geometry output.
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
  : public LocatedSetInterface
{ 
 public: 
  Polytope() { }
  Polytope(uint d) { }
  Polytope(const std::vector<Point>& v) : _vertices(v) { }

  void new_vertex(const Point& v) { 
    ARIADNE_ASSERT(this->_vertices.size()==0 || v.dimension()==this->_vertices.front().dimension()); 
    this->_vertices.push_back(v); }
  size_t number_of_vertices() const { return this->_vertices.size(); }
  const Point& vertex(size_t i) const { return this->_vertices[i]; }
  const std::vector<Point>& vertices() const { return this->_vertices; }
  std::vector<Point>& vertices() { return this->_vertices; }
  std::vector<Point>::const_iterator vertices_begin() const { return this->_vertices.begin(); }
  std::vector<Point>::const_iterator vertices_end() const { return this->_vertices.end(); }
 
  virtual Polytope* clone() const { return new Polytope(*this); }
  virtual uint dimension() const { if(this->_vertices.size()==0) { return 0; } return this->_vertices.front().dimension(); }
  virtual tribool disjoint(const Vector<Interval>& bx) const;
  virtual tribool intersects(const Vector<Interval>& bx) const;
  virtual tribool subset(const Vector<Interval>& bx) const;
  virtual Vector<Interval> bounding_box() const;
  virtual std::ostream& write(std::ostream& os) const { return os << *this; }
  
  friend Polytope convex_hull(const Polytope& p1, const Polytope& p2);
  friend Point baricentre(const Polytope& p);
 private:
  std::vector<Point> _vertices;
};

std::ostream& operator<<(std::ostream& os, const Polytope& p);




} // namespace Ariadne


#endif /* ARIADNE_POLYTOPE_H */
