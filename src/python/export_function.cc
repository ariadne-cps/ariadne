/***************************************************************************
 *            python/export_function.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/interpreted_function.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Python;

template<class R>
void export_function() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

  class_< InterpretedFunction<R> >("InterpretedFunction",init<>())
    .def(init< std::string >())
    .def("argument_size", &InterpretedFunction<R>::argument_size)
    .def("result_size", &InterpretedFunction<R>::result_size)
    .def("smoothness", &InterpretedFunction<R>::smoothness)
    .def("__call__",(Vector<A>(InterpretedFunction<R>::*)(const Vector<A>&)const)(&InterpretedFunction<R>::image))
    .def("jacobian",(Matrix<A>(InterpretedFunction<R>::*)(const Vector<A>&)const)(&InterpretedFunction<R>::jacobian))
    .def("read",(void(InterpretedFunction<R>::*)(const std::string&))(&InterpretedFunction<R>::read))
    .def(self_ns::str(self))
  ;
}

template void export_function<Rational>();
template void export_function<Float>();
