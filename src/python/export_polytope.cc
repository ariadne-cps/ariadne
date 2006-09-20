/***************************************************************************
 *            python/export_polytope.cc
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
#include "geometry/polytope.h"

#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;



void export_polytope() {
  typedef bool (*PolytPolytBinPred) (const RPolytope&, const RPolytope&);
  typedef RPolytope (*PolytPolytBinFunc) (const RPolytope&, const RPolytope&);

  def("disjoint", PolytPolytBinPred(&disjoint));
  def("interiors_intersect", PolytPolytBinPred(&interiors_intersect));
  def("inner_subset", PolytPolytBinPred(&inner_subset));
  def("subset", PolytPolytBinPred(&subset));
  def("convex_hull", PolytPolytBinFunc(&convex_hull));

  class_<RPolytope>("Polytope",init<size_type>())
    .def(init<RPointList>())
    .def(init<RPolytope>())
    .def(init<RRectangle>())
    .def("dimension", &RPolytope::dimension)
    .def("vertices", &RPolytope::vertices)
    .def("bounding_box", &RPolytope::bounding_box)
    .def(self_ns::str(self))
  ;

  
}
