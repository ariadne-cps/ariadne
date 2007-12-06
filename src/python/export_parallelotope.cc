/***************************************************************************
 *            python/export_parallelotope.cc
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

#include "python/float.h"

#include "geometry/parallelotope.h"

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"


#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class BS1, class BS2>
inline
BS1
touching_intersection(const BS1 &a, 
                      const BS2 &b) 
{
  if (disjoint(a,b)) {
    return BS1(a.dimension());
  }
  return a;
}

template<class R>
void export_parallelotope() 
{
  typedef Parallelotope<R> RParallelotope;
  typedef Zonotope<R> RZonotope;
  typedef Rectangle<R> RRectangle;
  typedef Point<R> RPoint;
  
  typedef LinearAlgebra::Matrix<R> RMatrix;
  
  def("touching_intersection",(RParallelotope(*)(const RParallelotope&,const RRectangle&))(&touching_intersection));
  def("touching_intersection", (RParallelotope(*)(const RParallelotope&,const RParallelotope&))(&touching_intersection));

  class_< RParallelotope, bases<RZonotope> > parallelotope_class("Parallelotope",init<int>());
  parallelotope_class.def(init<RPoint,RMatrix>());
  parallelotope_class.def(init<RParallelotope>());
  parallelotope_class.def(init<RRectangle>());
  parallelotope_class.def("empty", &RParallelotope::empty);
  parallelotope_class.def("contains", (tribool(RZonotope::*)(const RPoint&)const) &RZonotope::contains);
  parallelotope_class.def("divide", &RParallelotope::divide);
  parallelotope_class.def("subdivide", &RParallelotope::subdivide);
  parallelotope_class.def(self_ns::str(self));
  
/*
  // Can't use Python inheritence as the integration wrappers find Zonotope routine.
  class_<RParallelotope>("Parallelotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def("centre",&RZonotope::centre)
    .def("generators",&RZonotope::generators, return_value_policy<copy_const_reference>())
    .def("dimension", &RZonotope::dimension)
    .def("empty", &RZonotope::empty)
    .def("empty_interior", &RParallelotope::empty_interior)
    .def("contains", &RParallelotope::contains)
    .def("interior_contains", &RParallelotope::interior_contains)
    .def("divide", &RParallelotope::divide)
    .def("subdivide", &RParallelotope::subdivide)
    .def("bounding_box", &RZonotope::bounding_box)
    .def(self_ns::str(self))
  ;
*/
}

template void export_parallelotope<FloatPy>();
