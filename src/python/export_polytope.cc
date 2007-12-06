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

#include "python/float.h"

#include "geometry/rectangle.h"
#include "geometry/polytope.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_polytope() 
{
  typedef Polytope<R> RPolytope;

  def("disjoint", (tribool(*)(const RPolytope&, const RPolytope&))(&disjoint));
  def("subset", (tribool(*)(const RPolytope&, const RPolytope&))(&subset));
  def("convex_hull", (RPolytope(*)(const RPolytope&, const RPolytope&))(&convex_hull));

  class_< Polytope<R> >("Polytope",init<int>())
    .def(init< PointList<R> >())
    .def(init< Polytope<R> >())
    .def(init< Rectangle<R> >())
    .def("dimension", &Polytope<R>::dimension)
    .def("vertices", &Polytope<R>::vertices)
    .def("bounding_box", &Polytope<R>::bounding_box)
    .def(self_ns::str(self))
  ;

}

template void export_polytope<FloatPy>();
