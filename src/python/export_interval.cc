/***************************************************************************
 *            python/export_interval.cc
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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 

#include "numeric/interval.h"


#include <boost/python.hpp>
#include "python/python_utilities.h"
#include "python/python_float.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;


template<class R>
void export_interval()
{
  typedef Interval<R> I;

  class_< Interval<R> >(python_name<R>("Interval").c_str())
    //.def(init<std::string>())
    .def(init<int,int>())
    .def(init<double,double>())
    .def(init<R,R>())
    .def(init<int>())
    .def(init<double>())
    .def(init<R>())
    .def("__neg__", &neg<I,I>)
    .def("__add__", &add<I,I,int,I,R>)
    .def("__add__", &add<I,I,double,I,R>)
    .def("__add__", &add<I,I,R>)
    .def("__add__", &add<I,I,I>)
    .def("__radd__", &add<I,I,int>)
    .def("__radd__", &add<I,I,double>)
    .def("__radd__", &add<I,I,R>)
    .def("__sub__", &sub<I,I,int,I,R>)
    .def("__sub__", &sub<I,I,double,I,R>)
    .def("__sub__", &sub<I,I,R>)
    .def("__sub__", &sub<I,I,I>)
    .def("__rsub__", &rsub<I,I,int,I,R>)
    .def("__rsub__", &rsub<I,I,double,I,R>)
    .def("__rsub__", &rsub<I,I,R>)
    .def("__mul__", &mul<I,I,int,I,R>)
    .def("__mul__", &mul<I,I,double,I,R>)
    .def("__mul__", &mul<I,I,R>)
    .def("__mul__", &mul<I,I,I>)
    .def("__rmul__", &mul<I,I,int,I,R>)
    .def("__rmul__", &mul<I,I,double,I,R>)
    .def("__rmul__", &mul<I,I,R>)
    .def("__div__", &div<I,I,int>)
    .def("__div__", &div<I,I,double>)
    .def("__div__", &div<I,I,R>)
    .def("__div__", &div<I,I,I>)
    .def("__rdiv__", &rdiv<I,I,int,I,R>)
    .def("__rdiv__", &rdiv<I,I,double,I,R>)
    .def("__rdiv__", &rdiv<I,I,R>)
    .def("lower", &Interval<R>::lower, return_value_policy<copy_const_reference>())
    .def("upper", &Interval<R>::upper, return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
/*
  IFUN iabs(&boost::numeric::abs);
  IFUN iexp(&boost::numeric::exp);
  IFUN ilog(&boost::numeric::log);
  IFUN isin(&boost::numeric::sin);
  IFUN icos(&boost::numeric::cos);
  IFUN itan(&boost::numeric::tan);
  IFUN iasin(&boost::numeric::asin);
  IFUN iacos(&boost::numeric::acos);
  IFUN iatan(&boost::numeric::atan);
  IFUN isinh(&boost::numeric::sinh);
  IFUN icosh(&boost::numeric::cosh);
  IFUN itanh(&boost::numeric::tanh);
  IFUN iasinh(&boost::numeric::asinh);
  IFUN iacosh(&boost::numeric::acosh);
  IFUN iatanh(&boost::numeric::atanh);

  def("abs", iabs, "interval absolute value function");
  def("exp", iexp);
  def("log", ilog);
  def("sin", isin);
  def("cos", icos);
  def("tan", itan);
  def("asin", iasin);
  def("acos", iacos);
  def("atan", iatan);
  def("sinh", isinh);
  def("cosh", icosh);
  def("tanh", itanh);
  def("asinh", iasinh);
  def("acosh", iacosh);
  def("atanh", iatanh);
*/

}

template void export_interval<Float>();
