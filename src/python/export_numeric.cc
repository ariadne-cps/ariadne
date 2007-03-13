/***************************************************************************
 *            python/export_numeric.cc
 *
 *  22 June 2005
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


#include "python/python_utilities.h"

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"

#include "numeric/interval.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;

#include <boost/python.hpp>
using namespace boost::python;

typedef Interval<Float64> IFloat64;
typedef Interval<MPFloat> IMPFloat;

void export_numeric() {
  class_<Integer>("Integer")
    .def(init<int>())
    .def(init<Integer>())
    .def("__neg__", &neg<Integer,Integer>)
    .def("__add__", &add<Integer,Integer,Integer>)
    .def("__add__", &add<Integer,Integer,int>)
    .def("__radd__", &add<Integer,Integer,int>)
    .def("__sub__", &sub<Integer,Integer,Integer>)
    .def("__sub__", &sub<Integer,Integer,int>)
    .def("__rsub__", &rsub<Integer,Integer,int>)
    .def("__mul__", &mul<Integer,Integer,Integer>)
    .def("__mul__", &mul<Integer,Integer,int>)
    .def("__rmul__", &mul<Integer,Integer,int>)
    .def("__eq__", &eq<bool,Integer,int>)
    .def("__eq__", &eq<bool,Integer,Integer>)
    .def("__ne__", &ne<bool,Integer,int>)
    .def("__ne__", &ne<bool,Integer,Integer>)
    .def("__lt__", &lt<bool,Integer,int>)
    .def("__lt__", &lt<bool,Integer,Integer>)
    .def("__gt__", &gt<bool,Integer,int>)
    .def("__gt__", &gt<bool,Integer,Integer>)
    .def("__le__", &le<bool,Integer,int>)
    .def("__le__", &le<bool,Integer,Integer>)
    .def("__ge__", &ge<bool,Integer,int>)
    .def("__ge__", &ge<bool,Integer,Integer>)
    .def(self_ns::str(self))
  ;

  class_<Float64>("Float64")
    .def(init<int>())
    .def(init<double>())
    .def(init<Float64>())
    .def("__neg__", &neg<Float64,Float64>)
    .def("__add__", &add<IFloat64,Float64,int,Float64,Float64>)
    .def("__add__", &add<IFloat64,Float64,double,Float64,Float64>)
    .def("__add__", &add<IFloat64,Float64,Float64>)
    .def("__radd__", &add<IFloat64,Float64,int,Float64,Float64>)
    .def("__radd__", &add<IFloat64,Float64,double,Float64,Float64>)
    .def("__sub__", &sub<IFloat64,Float64,int,Float64,Float64>)
    .def("__sub__", &sub<IFloat64,Float64,double,Float64,Float64>)
    .def("__sub__", &sub<IFloat64,Float64,Float64>)
    .def("__rsub__", &rsub<IFloat64,Float64,int,Float64,Float64>)
    .def("__rsub__", &rsub<IFloat64,Float64,double,Float64,Float64>)
    .def("__mul__", &mul<IFloat64,Float64,int,Float64,Float64>)
    .def("__mul__", &mul<IFloat64,Float64,double,Float64,Float64>)
    .def("__mul__", &mul<IFloat64,Float64,Float64>)
    .def("__rmul__", &mul<IFloat64,Float64,int,Float64,Float64>)
    .def("__rmul__", &mul<IFloat64,Float64,double,Float64,Float64>)
    .def("__div__", &div<IFloat64,Float64,int,Float64,Float64>)
    .def("__div__", &div<IFloat64,Float64,double,Float64,Float64>)
    .def("__div__", &div<IFloat64,Float64,Float64>)
    .def("__rdiv__", &rdiv<IFloat64,Float64,int,Float64,Float64>)
    .def("__rdiv__", &rdiv<IFloat64,Float64,double,Float64,Float64>)
    .def("__eq__", &eq<bool,Float64,double>)
    .def("__eq__", &eq<bool,Float64,Float64>)
    .def("__ne__", &ne<bool,Float64,double>)
    .def("__ne__", &ne<bool,Float64,Float64>)
    .def("__lt__", &lt<bool,Float64,double>)
    .def("__lt__", &lt<bool,Float64,Float64>)
    .def("__gt__", &gt<bool,Float64,double>)
    .def("__gt__", &gt<bool,Float64,Float64>)
    .def("__le__", &le<bool,Float64,double>)
    .def("__le__", &le<bool,Float64,Float64>)
    .def("__ge__", &ge<bool,Float64,double>)
    .def("__ge__", &ge<bool,Float64,Float64>)
    .def(self_ns::str(self))
  ;

  class_<MPFloat>("MPFloat")
    .def(init<int>())
    .def(init<Integer>())
    .def(init<double>())
    .def(init<MPFloat>())
    .def("__neg__", &neg<MPFloat,MPFloat>)
    .def("__add__", &add<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__add__", &add<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__add__", &add<IMPFloat,MPFloat,MPFloat>)
    .def("__radd__", &add<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__radd__", &add<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__sub__", &sub<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__sub__", &sub<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__sub__", &sub<IMPFloat,MPFloat,MPFloat>)
    .def("__rsub__", &rsub<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__rsub__", &rsub<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__mul__", &mul<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__mul__", &mul<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__mul__", &mul<IMPFloat,MPFloat,MPFloat>)
    .def("__rmul__", &mul<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__rmul__", &mul<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__div__", &div<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__div__", &div<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__div__", &div<IMPFloat,MPFloat,MPFloat>)
    .def("__rdiv__", &rdiv<IMPFloat,MPFloat,int,MPFloat,MPFloat>)
    .def("__rdiv__", &rdiv<IMPFloat,MPFloat,double,MPFloat,MPFloat>)
    .def("__eq__", &eq<bool,MPFloat,double>)
    .def("__eq__", &eq<bool,MPFloat,MPFloat>)
    .def("__ne__", &ne<bool,MPFloat,double>)
    .def("__ne__", &ne<bool,MPFloat,MPFloat>)
    .def("__lt__", &lt<bool,MPFloat,double>)
    .def("__lt__", &lt<bool,MPFloat,MPFloat>)
    .def("__gt__", &gt<bool,MPFloat,double>)
    .def("__gt__", &gt<bool,MPFloat,MPFloat>)
    .def("__le__", &le<bool,MPFloat,double>)
    .def("__le__", &le<bool,MPFloat,MPFloat>)
    .def("__ge__", &ge<bool,MPFloat,double>)
    .def("__ge__", &ge<bool,MPFloat,MPFloat>)
    .def("precision", &MPFloat::precision)
    .def(self_ns::str(self))
  ;

  class_<Rational>("Rational")
    .def(init<int,int>())
    .def(init<Integer,Integer>())
    .def(init<int>())
    .def(init<Integer>())
    .def(init<double>())
    .def(init<Float64>())
    .def(init<MPFloat>())
    .def(init<Rational>())
    .def("__neg__", &neg<Rational,Rational>)
    .def("__add__", &add<Rational,Rational,int>)
    .def("__add__", &add<Rational,Rational,Integer>)
    .def("__add__", &add<Rational,Rational,double>)
    .def("__add__", &add<Rational,Rational,Float64,Rational,Rational>)
    .def("__add__", &add<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__add__", &add<Rational,Rational,Rational>)
    .def("__radd__", &add<Rational,Rational,int>)
    .def("__radd__", &add<Rational,Rational,Integer>)
    .def("__radd__", &add<Rational,Rational,double>)
    .def("__radd__", &add<Rational,Rational,Float64,Rational,Rational>)
    .def("__radd__", &add<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__sub__", &sub<Rational,Rational,int>)
    .def("__sub__", &sub<Rational,Rational,Integer>)
    .def("__sub__", &sub<Rational,Rational,double>)
    .def("__sub__", &sub<Rational,Rational,Float64,Rational,Rational>)
    .def("__sub__", &sub<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__sub__", &sub<Rational,Rational,Rational>)
    .def("__rsub__", &rsub<Rational,Rational,int>)
    .def("__rsub__", &rsub<Rational,Rational,Integer>)
    .def("__rsub__", &rsub<Rational,Rational,double>)
    .def("__rsub__", &rsub<Rational,Rational,Float64,Rational,Rational>)
    .def("__rsub__", &rsub<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__mul__", &mul<Rational,Rational,int>)
    .def("__mul__", &mul<Rational,Rational,Integer>)
    .def("__mul__", &mul<Rational,Rational,double>)
    .def("__mul__", &mul<Rational,Rational,Float64,Rational,Rational>)
    .def("__mul__", &mul<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__mul__", &mul<Rational,Rational,Rational>)
    .def("__rmul__", &mul<Rational,Rational,int>)
    .def("__rmul__", &mul<Rational,Rational,Integer>)
    .def("__rmul__", &mul<Rational,Rational,double>)
    .def("__rmul__", &mul<Rational,Rational,Float64,Rational,Rational>)
    .def("__rmul__", &mul<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__div__", &div<Rational,Rational,int>)
    .def("__div__", &div<Rational,Rational,Integer>)
    .def("__div__", &div<Rational,Rational,double>)
    .def("__div__", &div<Rational,Rational,Float64,Rational,Rational>)
    .def("__div__", &div<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__div__", &div<Rational,Rational,Rational>)
    .def("__rdiv__", &rdiv<Rational,Rational,int>)
    .def("__rdiv__", &rdiv<Rational,Rational,Integer>)
    .def("__rdiv__", &rdiv<Rational,Rational,double>)
    .def("__rdiv__", &rdiv<Rational,Rational,Float64,Rational,Rational>)
    .def("__rdiv__", &rdiv<Rational,Rational,MPFloat,Rational,Rational>)
    .def("__eq__", &eq<bool,Rational,Rational>)
    .def("__eq__", &eq<bool,Rational,double>)
    .def("__eq__", &eq<bool,Rational,Integer>)
    .def("__eq__", &eq<bool,Rational,double>)
    .def("__ne__", &ne<bool,Rational,Rational>)
    .def("__ne__", &ne<bool,Rational,double>)
    .def("__ne__", &ne<bool,Rational,Integer>)
    .def("__ne__", &ne<bool,Rational,double>)
    .def("__lt__", &lt<bool,Rational,Rational>)
    .def("__lt__", &lt<bool,Rational,double>)
    .def("__lt__", &lt<bool,Rational,Integer>)
    .def("__lt__", &lt<bool,Rational,double>)
    .def("__gt__", &gt<bool,Rational,double>)
    .def("__gt__", &gt<bool,Rational,Integer>)
    .def("__gt__", &gt<bool,Rational,double>)
    .def("__le__", &le<bool,Rational,Rational>)
    .def("__le__", &le<bool,Rational,double>)
    .def("__le__", &le<bool,Rational,Integer>)
    .def("__le__", &le<bool,Rational,double>)
    .def("__ge__", &ge<bool,Rational,double>)
    .def("__ge__", &ge<bool,Rational,Integer>)
    .def("__ge__", &ge<bool,Rational,double>)
    .def("numerator", &Rational::numerator)
    .def("denominator", &Rational::denominator)
    .def(self_ns::str(self))
  ;
 
  
}
