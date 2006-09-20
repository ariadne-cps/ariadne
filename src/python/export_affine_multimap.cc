/***************************************************************************
 *            python/export_affine_map_with_set.cc
 *
 *  31 May 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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


#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/simplex.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "system/affine_multimap.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

typedef RRectangle (RAffineMultiMapRect::*AffMultiMapRectCallPoint) (const RPoint&) const;
typedef RRectangle (RAffineMultiMapRect::*AffMultiMapRectCallRectangle) (const RRectangle&) const;
//typedef RParallelotope (RAffineMultiMapRect::*AffMultiMapRectCallParallelotope) (const RParallelotope&) const;
//typedef RZonotope (RAffineMultiMapRect::*AffMultiMapRectCallZonotope) (const RZonotope&) const;
//typedef RSimplex (RAffineMultiMapRect::*AffMultiMapRectCallSimplex) (const RSimplex&) const;
//typedef RPolytope (RAffineMultiMapRect::*AffMultiMapRectCallPolytope) (const RPolytope&) const;
//typedef RParallelotopeListSet (RAffineMultiMapRect::*AffMultiMapRectCallParallelotopeListSet) (const RParallelotopeListSet&) const;
//typedef RZonotopeListSet (RAffineMultiMapRect::*AffMultiMapRectCallZonotopeListSet) (const RZonotopeListSet&) const;
//typedef RZonotopeListSet (RAffineMultiMapRect::*AffMultiMapRectCallGridMaskSet) (const RGridMaskSet&) const;

AffMultiMapRectCallPoint affine_multimap_rect_call_point=&RAffineMultiMapRect::operator();
AffMultiMapRectCallRectangle affine_multimap_rect_call_rectangle=&RAffineMultiMapRect::operator();
//AffMultiMapRectCallParallelotope affine_multimap_rect_call_parallelotope=&RAffineMultiMapRect::operator();
//AffMultiMapRectCallSimplex affine_multimap_rect_call_simplex=&RAffineMultiMapRect::operator();
//AffMultiMapRectCallZonotope affine_multimap_rect_call_zonotope=&RAffineMultiMapRect::operator();
//AffMultiMapRectCallPolytope affine_multimap_rect_call_polytope=&RAffineMultiMapRect::operator();
//AffMultiMapRectCallParallelotopeListSet affine_multimap_rect_call_parallelotope_list_set=&RAffineMultiMapRect::operator();
//AffMultiMapRectCallZonotopeListSet affine_multimap_rect_call_zonotope_list_set=&RAffineMultiMapRect::operator();
//AffMultiMapRectCallGridMaskSet affine_multimap_rect_call_grid_mask_set=&RAffineMultiMapRect::operator();

//typedef RParallelotope (RAffineMultiMapPltp::*AffMultiMapPltpCallPoint) (const RPoint&) const;
//typedef RParallelotope (RAffineMultiMapPltp::*AffMultiMapPltpCallParallelotope) (const RParallelotope&) const;
//typedef RZonotope (RAffineMultiMapPltp::*AffMultiMapPltpCallZonotope) (const RZonotope&) const;
//typedef RPolytope (RAffineMultiMapPltp::*AffMultiMapPltpCallPolytope) (const RPolytope&) const;
//typedef RParallelotopeListSet (RAffineMultiMapPltp::*AffMultiMapPltpCallParallelotopeListSet) (const RParallelotopeListSet&) const;
//typedef RZonotopeListSet (RAffineMultiMapPltp::*AffMultiMapPltpCallZonotopeListSet) (const RZonotopeListSet&) const;
//typedef RZonotopeListSet (RAffineMultiMapPltp::*AffMultiMapPltpCallGridMaskSet) (const RGridMaskSet&) const;

//AffMultiMapPltpCallPoint affine_multimap_pltp_call_point=&RAffineMultiMapPltp::operator();
//AffMultiMapPltpCallParallelotope affine_multimap_pltp_call_parallelotope=&RAffineMultiMapPltp::operator();
//AffMultiMapPltpCallZonotope affine_multimap_pltp_call_zonotope=&RAffineMultiMapPltp::operator();
//AffMultiMapPltpCallPolytope affine_multimap_pltp_call_polytope=&RAffineMultiMapPltp::operator();
//AffMultiMapPltpCallParallelotopeListSet affine_multimap_pltp_call_parallelotope_list_set=&RAffineMultiMapPltp::operator();
//AffMultiMapPltpCallZonotopeListSet affine_multimap_pltp_call_zonotope_list_set=&RAffineMultiMapPltp::operator();
//AffMultiMapPltpCallGridMaskSet affine_multimap_pltp_call_grid_mask_set=&RAffineMultiMapPltp::operator();

typedef RZonotope (RAffineMultiMapZntp::*AffMultiMapZntpCallPoint) (const RPoint&) const;
typedef RZonotope (RAffineMultiMapZntp::*AffMultiMapZntpCallZonotope) (const RZonotope&) const;
//typedef RPolytope (RAffineMultiMapZntp::*AffMultiMapZntpCallPolytope) (const RPolytope&) const;
//typedef RParallelotopeListSet (RAffineMultiMapZntp::*AffMultiMapZntpCallParallelotopeListSet) (const RParallelotopeListSet&) const;
//typedef RZonotopeListSet (RAffineMultiMapZntp::*AffMultiMapZntpCallZonotopeListSet) (const RZonotopeListSet&) const;
//typedef RZonotopeListSet (RAffineMultiMapZntp::*AffMultiMapZntpCallGridMaskSet) (const RGridMaskSet&) const;

AffMultiMapZntpCallPoint affine_multimap_zntp_call_point=&RAffineMultiMapZntp::operator();
AffMultiMapZntpCallZonotope affine_multimap_zntp_call_zonotope=&RAffineMultiMapZntp::operator();
//AffMultiMapZntpCallPolytope affine_multimap_zntp_call_polytope=&RAffineMultiMapZntp::operator();
//AffMultiMapZntpCallZonotopeListSet affine_multimap_zntp_call_zonotope_list_set=&RAffineMultiMapZntp::operator();
//AffMultiMapZntpCallGridMaskSet affine_multimap_zntp_call_grid_mask_set=&RAffineMultiMapZntp::operator();

void export_affine_map_with_set() {

  class_< RAffineMultiMapRect>("AffineMultiMap",init<RMatrix,RRectangle>())
    .def(init<RAffineMultiMapRect>())
    .def("argument_dimension", &RAffineMultiMapRect::argument_dimension)
    .def("result_dimension", &RAffineMultiMapRect::result_dimension)
    .def("__call__", affine_multimap_rect_call_point)
    .def("__call__", affine_multimap_rect_call_rectangle)
//    .def("__call__", affine_multimap_rect_call_polytope)
//    .def("__call__", affine_multimap_rect_call_zonotope)
//    .def("__call__", affine_multimap_rect_call_zonotope_list_set)
//    .def("__call__", affine_multimap_rect_call_grid_mask_set)
//    .def(self_ns::str(self))    // __self_ns::str__
  ;

/*
  class_< RAffineMultiMapPltp, bases<RAffineMap> >("AffineMultiMapParallelotope",init<RMatrix,RParallelotope>())
    .def(init<RAffineMultiMapPltp>())
    .def("argument_dimension", &RAffineMultiMapPltp::argument_dimension)
    .def("result_dimension", &RAffineMultiMapPltp::result_dimension)
    .def("__call__", affine_multimap_pltp_call_point)
    .def("__call__", affine_multimap_pltp_call_parallelotope)
    .def("__call__", affine_multimap_pltp_call_polytope)
    .def("__call__", affine_multimap_pltp_call_zonotope)
    .def("__call__", affine_multimap_pltp_call_parallelotope_list_set)
    .def("__call__", affine_multimap_pltp_call_zonotope_list_set)
    .def("__call__", affine_multimap_pltp_call_grid_mask_set)
//    .def(self_ns::str(self))    // __self_ns::str__
  ;
*/

  class_< RAffineMultiMapZntp>("AffineMultiMapZonotope",init<RMatrix,RZonotope>())
    .def(init<RAffineMultiMapZntp>())
    .def("argument_dimension", &RAffineMultiMapZntp::argument_dimension)
    .def("result_dimension", &RAffineMultiMapZntp::result_dimension)
    .def("__call__", affine_multimap_zntp_call_point)
//    .def("__call__", affine_multimap_zntp_call_polytope)
    .def("__call__", affine_multimap_zntp_call_zonotope)
//    .def("__call__", affine_multimap_zntp_call_zonotope_list_set)
//    .def("__call__", affine_multimap_zntp_call_grid_mask_set)
//    .def(self_ns::str(self))    // __self_ns::str__
  ;

}
