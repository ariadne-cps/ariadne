/***************************************************************************
 *            python/export_polyhedron.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <vector>


#include "utility/stlio.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/polyhedron.h"

#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;


typedef std::vector<RPoint> RPointList;

inline RPoint point_list_get(const RPointList& pl, const size_type& n) {
  return pl[n];
}

void export_polyhedron() {
  typedef bool (*PolyPolyBinPred) (const RPolyhedron&, const RPolyhedron&);
  typedef RPolyhedron (*PolyPolyBinFunc) (const RPolyhedron&, const RPolyhedron&);

  def("intersection", PolyPolyBinFunc(&intersection));
  def("regular_intersection", PolyPolyBinFunc(&regular_intersection));
  def("convex_hull", PolyPolyBinFunc(&convex_hull));
  def("interiors_intersect", PolyPolyBinPred(&interiors_intersect));
  def("disjoint", PolyPolyBinPred(&disjoint));
  def("inner_subset", PolyPolyBinPred(&inner_subset));
  def("subset", PolyPolyBinPred(&subset));

  class_<RPolyhedron>("Polyhedron",init<int>())
    .def(init<RMatrix,RVector>())
    .def(init<RPointList>())
    .def(init<RPolyhedron>())
    .def(init<RRectangle>())
    .def(init<RParallelotope>())
    .def("dimension", &RPolyhedron::dimension)
    .def("vertices", &RPolyhedron::vertices)
    .def("bounding_box", &RPolyhedron::bounding_box)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  typedef RPoint (RPointList::* RPLGetFunc) (const RPointList::size_type&) const;
  class_<RPointList>("PointList",init<>())
    .def("size", &RPointList::size)
    .def("append", &RPointList::push_back)
    .def("__getitem__", &point_list_get)
   // .def(self_ns::str(self))    // __self_ns::str__
  ;

  
}
