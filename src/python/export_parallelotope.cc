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



#include "geometry/parallelotope.h"

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template <template <typename> class BS1, template<typename> class BS2>
inline
BS1<Real> 
touching_intersection(const BS1<Real> &a, 
                      const BS2<Real> &b) 
{
  if (disjoint(a,b))
    return BS1<Real>(a.dimension());
  
  return a;
}

template Parallelotope<Real> touching_intersection(
                const Parallelotope<Real> &,  
                const Rectangle<Real> &);

template Parallelotope<Real> touching_intersection(
                const Parallelotope<Real> &,  
                const Parallelotope<Real> &);

void export_parallelotope() {
  typedef RParallelotope (*PltpRectBinFun) (const RParallelotope&, 
                                            const RRectangle&);
  typedef RParallelotope (*PltpPltpBinFun) (const RParallelotope&, 
                                            const RParallelotope&);
  typedef bool (*PltpPltpBinPred) (const RParallelotope&, 
                                   const RParallelotope&);
  typedef bool (*ZntpZntpBinPred) (const RZonotope&, 
                                   const RZonotope&);

  def("touching_intersection", PltpRectBinFun(&touching_intersection));
  def("touching_intersection", PltpPltpBinFun(&touching_intersection));

/*
  class_< RParallelotope, bases<RZonotope> >("Parallelotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("empty", &RParallelotope::empty)
    .def("empty_interior", &RParallelotope::empty_interior)
    .def("contains", &RParallelotope::contains)
    .def("interior_contains", &RParallelotope::interior_contains)
    .def("divide", &RParallelotope::divide)
    .def("subdivide", &RParallelotope::subdivide)
    .def(self_ns::str(self))
  ;
*/
  // Can't use Python inheritence as the integration wrappers find Zonotope routine.
  class_<RParallelotope>("Parallelotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
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
  
}
