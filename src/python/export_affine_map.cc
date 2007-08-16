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

#include "python/python_float.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/simplex.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "system/affine_map.h"


#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_affine_map() 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  typedef Interval<R> I;
  
  class_< AffineMap<R>, bases< MapInterface<R> > >("AffineMap",init< Matrix<R>, Vector<R> >())
    .def("argument_dimension",&AffineMap<R>::argument_dimension)
    .def("result_dimension",&AffineMap<R>::result_dimension)
    .def("smoothness",&AffineMap<R>::smoothness)
    .def("__call__",(Point<F>(AffineMap<R>::*)(const Point<F>&)const)(&AffineMap<R>::image))
    .def("__call__",(Zonotope<F>(AffineMap<R>::*)(const Zonotope<F>&)const)(&AffineMap<R>::image))
    .def("__call__",(Polytope<F>(AffineMap<R>::*)(const Polytope<F>&)const)(&AffineMap<R>::image))
    .def(self_ns::str(self))
  ;
}

template void export_affine_map<Float>();
