/***************************************************************************
 *            python/export_map.cc
 *
 *  13 February 2006
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

#include "python/python_float.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "system/function.h"
#include "system/map.h"
#include "system/function_map.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
class MapWrapper : public MapInterface<R>, public wrapper< MapInterface<R> >
{
  typedef typename MapInterface<R>::F F;
 public:
  MapInterface<R>* clone() const { return this->get_override("clone")(); }
  Point<F> image(const Point<F>&) const { return this->get_override("clone")(); }
  dimension_type argument_dimension() const { return this->get_override("argument_dimension")(); }
  dimension_type result_dimension() const { return this->get_override("result_dimension")(); }
  size_type smoothness() const { return this->get_override("smoothness")(); }
  std::string name() const { return this->get_override("name")(); }
};

template<class R>
void export_map() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

  class_<MapWrapper<R>, boost::noncopyable>("MapInterface")
    .def("argument_dimension", pure_virtual(&MapWrapper<R>::argument_dimension))
    .def("result_dimension", pure_virtual(&MapWrapper<R>::result_dimension))
    .def("smoothness", pure_virtual(&MapWrapper<R>::smoothness))
  ;

  class_< FunctionMap<R>, bases< MapInterface<R> > >("FunctionMap", init< const Function<R>&, Point<A> >())
    .def("argument_dimension", &FunctionMap<R>::argument_dimension)
    .def("result_dimension", &FunctionMap<R>::result_dimension)
    .def("number_of_parameters", &FunctionMap<R>::number_of_parameters)
    .def("smoothness", &FunctionMap<R>::smoothness)
    .def("parameters", &FunctionMap<R>::parameters, return_value_policy<copy_const_reference>())
    .def("set_parameters", &FunctionMap<R>::set_parameters)
    .def("__call__",(Point<A>(FunctionMap<R>::*)(const Point<A>&)const)(&FunctionMap<R>::image))
    .def("image",(Point<A>(FunctionMap<R>::*)(const Point<A>&)const)(&FunctionMap<R>::image))
    .def("jacobian",(Matrix<A>(FunctionMap<R>::*)(const Point<A>&)const)(&FunctionMap<R>::jacobian))
    .def(self_ns::str(self))
  ;
}

template void export_map<Float>();
