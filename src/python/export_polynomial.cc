/***************************************************************************
 *            python/export_polynomial.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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
#include "function/polynomial.h"
#include "geometry/box.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_polynomial() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

  class_< Polynomial<R> >("Polynomial",init<>())
    .def(init< Polynomial<R> >())
    .def("result_size", &Polynomial<R>::result_size)
    .def("argument_size", &Polynomial<R>::argument_size)
    .def("degree", &Polynomial<R>::degree)
    .def("set",(void(Polynomial<R>::*)(const size_type&,const MultiIndex&,const R&)) &Polynomial<R>::set)
    .def("get",(R(Polynomial<R>::*)(const size_type&,const MultiIndex&)const) &Polynomial<R>::get)
    .def("component",(Polynomial<R>(Polynomial<R>::*)(const size_type&)const) &Polynomial<R>::component)
    .def("truncate",(Polynomial<R>(Polynomial<R>::*)(const size_type&,const size_type&,const Box<R>&)const) &Polynomial<R>::truncate)
    .def("evaluate",(Vector<R>(Polynomial<R>::*)(const Vector<R>&)const) &Polynomial<R>::evaluate)
    .def("jacobian",(Matrix<R>(Polynomial<R>::*)(const Vector<R>&)const) &Polynomial<R>::jacobian)
    .def("detivative",(Polynomial<R>(Polynomial<R>::*)(const size_type&)const) &Polynomial<R>::derivative)
    .def("__add__",(Polynomial<A>(*)(const Polynomial<R>&,const Polynomial<R>&)) &operator+)
    .def("__sub__",(Polynomial<A>(*)(const Polynomial<R>&,const Polynomial<R>&)) &operator-)
    .def("__mul__",(Polynomial<A>(*)(const Polynomial<R>&,const Polynomial<R>&)) &operator*)
    .def("pow",(Polynomial<A>(*)(const Polynomial<R>&,const uint&)) &pow)
    .def(self_ns::str(self))
  ;
 
}

template void export_polynomial<FloatPy>();
