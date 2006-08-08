/***************************************************************************
 *            python/export_affine_map.cc
 *
 *  6 February 2006
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
#include "geometry/polyhedron.h"
#include "system/affine_map.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

typedef System::Map<Real> RMap;

typedef RPoint (RAffineMap::*AffMapCallPoint) (const RPoint&) const;
typedef RParallelotope (RAffineMap::*AffMapCallParallelotope) (const RParallelotope&) const;
typedef RRectangle (RAffineMap::*AffMapCallRectangle) (const RRectangle&) const;
typedef RZonotope (RAffineMap::*AffMapCallZonotope) (const RZonotope&) const;
typedef RSimplex (RAffineMap::*AffMapCallSimplex) (const RSimplex&) const;
typedef RPolyhedron (RAffineMap::*AffMapCallPolyhedron) (const RPolyhedron&) const;
typedef RPolyhedron (RMap::*MapCallPolyhedron) (const RPolyhedron&) const;
typedef RSimplex (RMap::*MapCallSimplex) (const RSimplex&) const;

AffMapCallPoint affine_map_call_point=&RAffineMap::operator();
AffMapCallRectangle affine_map_call_rectangle=&RAffineMap::operator();
AffMapCallParallelotope affine_map_call_parallelotope=&RAffineMap::operator();
AffMapCallZonotope affine_map_call_zonotope=&RAffineMap::operator();
MapCallSimplex affine_map_call_simplex=&RMap::operator();
MapCallPolyhedron affine_map_call_polyhedron=&RMap::operator();


void export_affine_map() {

  class_< RAffineMap, bases<RMapBase> >("AffineMap",init<RMatrix,RVector>())
    .def("argument_dimension", &RAffineMap::argument_dimension)
    .def("result_dimension", &RAffineMap::result_dimension)
    .def("__call__", affine_map_call_point)
    .def("__call__", affine_map_call_rectangle)
    .def("__call__", affine_map_call_parallelotope)
    .def("__call__", affine_map_call_zonotope)
    .def("__call__", affine_map_call_simplex)
    .def("__call__", affine_map_call_polyhedron)
    .def(self_ns::str(self))
  ;
}
