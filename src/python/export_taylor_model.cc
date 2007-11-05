/***************************************************************************
 *            python/export_taylor_model.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

#include "linear_algebra/vector.h"
#include "function/taylor_model.h"
#include "geometry/rectangle.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_taylor_model() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

  class_< TaylorModel<R> >("TaylorModel",init<>())
    .def(init< TaylorModel<R> >())
    .def("result_size", &TaylorModel<R>::result_size)
    .def("argument_size", &TaylorModel<R>::argument_size)
    .def("order", &TaylorModel<R>::order)
    .def("smoothness", &TaylorModel<R>::smoothness)
    .def("set",(void(TaylorModel<R>::*)(const size_type&,const MultiIndex&,const R&)) &TaylorModel<R>::set)
    .def("get",(R(TaylorModel<R>::*)(const size_type&,const MultiIndex&)const) &TaylorModel<R>::get)
    .def("component",(TaylorModel<R>(TaylorModel<R>::*)(const size_type&)const) &TaylorModel<R>::component)
    .def("truncate",(TaylorModel<R>(TaylorModel<R>::*)(const size_type&,const size_type&,const Rectangle<R>&)const) &TaylorModel<R>::truncate)
    .def("evaluate",(Vector<R>(TaylorModel<R>::*)(const Vector<R>&)const) &TaylorModel<R>::evaluate)
    .def("jacobian",(Matrix<R>(TaylorModel<R>::*)(const Vector<R>&)const) &TaylorModel<R>::jacobian)
    .def("derivative",(TaylorModel<R>(TaylorModel<R>::*)(const size_type&)const) &TaylorModel<R>::derivative)
    .def("__add__",(TaylorModel<A>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &operator+)
    .def("__sub__",(TaylorModel<A>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &operator-)
    .def("__mul__",(TaylorModel<A>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &operator*)
    .def("__div__",(TaylorModel<A>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &operator*)
    .def("pow",(TaylorModel<A>(*)(const TaylorModel<R>&,const uint&)) &pow)
    .def(self_ns::str(self))
  ;
 
  def("compose",(TaylorModel<A>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &compose);
  def("inverse",(TaylorModel<A>(*)(const TaylorModel<R>&,const Vector<R>&)) &inverse);
  def("implicit",(TaylorModel<A>(*)(const TaylorModel<R>&,const Vector<R>&)) &implicit);
}

template void export_taylor_model<Float>();
