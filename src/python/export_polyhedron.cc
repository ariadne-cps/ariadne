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

#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"

#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;



void export_polyhedron() {
  typedef bool (*PolyhPolyhBinPred) (const RPolyhedron&, const RPolyhedron&);
  typedef RPolyhedron (*PolyhPolyhBinFunc) (const RPolyhedron&, const RPolyhedron&);

  def("disjoint", PolyhPolyhBinPred(&disjoint));
  def("interiors_intersect", PolyhPolyhBinPred(&interiors_intersect));
  def("inner_subset", PolyhPolyhBinPred(&inner_subset));
  def("subset", PolyhPolyhBinPred(&subset));
  def("convex_hull", PolyhPolyhBinFunc(&convex_hull));

  class_<RPolyhedron>("Polyhedron",init<int>())
    .def(init<RPointList>())
    .def(init<RPolyhedron>())
    .def(init<RRectangle>())
    .def("dimension", &RPolyhedron::dimension)
    .def("vertices", &RPolyhedron::vertices)
    .def("bounding_box", &RPolyhedron::bounding_box)
    .def(self_ns::str(self))
  ;
  

  
}
