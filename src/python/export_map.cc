/***************************************************************************
 *            python/export_map.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/taylor_derivative.h"
#include "function/function_interface.h"
#include "function/affine_function.h"
#include "geometry/point.h"
#include "geometry/box.h"

#include "system/map.h"
#include "system/affine_map.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_map() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;


  class_< Map<R> >("Map", init< const FunctionInterface<R>& >())
    .def("argument_dimension", &Map<R>::argument_dimension)
    .def("result_dimension", &Map<R>::result_dimension)
    .def("smoothness", &Map<R>::smoothness)
    //.def("number_of_parameters", &Map<R>::number_of_parameters)
    //.def("parameters", &Map<R>::parameters, return_value_policy<copy_const_reference>())
    .def("__call__",(Point<A>(Map<R>::*)(const Point<A>&)const)(&Map<R>::image))
    .def("image",(Point<A>(Map<R>::*)(const Point<A>&)const)(&Map<R>::image))
    .def("jacobian",(Matrix<A>(Map<R>::*)(const Point<A>&)const)(&Map<R>::jacobian))
    .def("derivative",(TaylorDerivative<A>(Map<R>::*)(const Point<A>&,const smoothness_type&)const)(&Map<R>::derivative))
    .def(self_ns::str(self))
  ;

  class_< AffineMap<R>, bases< Map<R> > >("AffineMap", init<const Matrix<R>&,const Vector<R>&>())
    .def(self_ns::str(self));


}

template void export_map<FloatPy>();
