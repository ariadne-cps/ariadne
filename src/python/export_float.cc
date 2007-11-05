/***************************************************************************
 *            python/export_float.cc
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


#include <iostream>
#include <iomanip>

#include "python/python_utilities.h"

#include "python/python_float.h"
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/interval.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

void set_output_precision(uint);
template<class R> void set_default_precision(uint);
template<class R> uint default_precision();

void set_output_precision(uint p) {
  std::cout << std::setprecision(p);
}

#if PYTHON_FLOAT == Float64 

template<> void set_default_precision<Float64>(uint p) { 
  throw std::runtime_error("Cannot set precision of 64-bit float");
}

template<> uint default_precision<Float64>() { 
  return 64;
}

#elif PYTHON_FLOAT == Float64 

template<> void set_default_precision<FloatMP>(uint p) { 
  FloatMP::set_default_precision(p);
}

template<> uint default_precision<FloatMP>() { 
  return FloatMP::default_precision();
}

#else

#endif

template<class R> 
std::string
__str__(const R& x) {
  std::stringstream ss;
  if(x==floor(x)) {
    ss << x << ".";
  } else {
    //ss << std::fixed;
    ss << std::setprecision(20);
    ss << x;
  }
  return ss.str();
}

template<class R> 
std::string
__repr__(const R& x) {
  std::stringstream ss;
  ss << "Float(";
  if(x==floor(x)) {
    ss << x;
  } else {
    //ss << std::fixed;
    ss << std::setprecision(20);
    ss << x;
  }
  ss << ")";
  return ss.str();
}

template<class R>
void 
export_float() 
{
  typedef R Float;
  typedef Interval<R> IFloat;
  
  def("set_output_precision",(void(*)(uint))&set_output_precision);

  class_<Float>("Float")
    //.def("__init__", make_constructor(&make_float<R>) )
    .def(init<int>())
    .def(init<double>())
    .def(init<std::string>())
    .def(init<Float>())
    .def("__neg__", &neg<Float,Float>)
    .def("__add__", &add<IFloat,Float,int,Float,Float>)
    .def("__add__", &add<IFloat,Float,double,Float,Float>)
    .def("__add__", &add<IFloat,Float,Float>)
    .def("__radd__", &add<IFloat,Float,int,Float,Float>)
    .def("__radd__", &add<IFloat,Float,double,Float,Float>)
    .def("__sub__", &sub<IFloat,Float,int,Float,Float>)
    .def("__sub__", &sub<IFloat,Float,double,Float,Float>)
    .def("__sub__", &sub<IFloat,Float,Float>)
    .def("__rsub__", &rsub<IFloat,Float,int,Float,Float>)
    .def("__rsub__", &rsub<IFloat,Float,double,Float,Float>)
    .def("__mul__", &mul<IFloat,Float,int,Float,Float>)
    .def("__mul__", &mul<IFloat,Float,double,Float,Float>)
    .def("__mul__", &mul<IFloat,Float,Float>)
    .def("__rmul__", &mul<IFloat,Float,int,Float,Float>)
    .def("__rmul__", &mul<IFloat,Float,double,Float,Float>)
    .def("__div__", &div<IFloat,Float,int,Float,Float>)
    .def("__div__", &div<IFloat,Float,double,Float,Float>)
    .def("__div__", &div<IFloat,Float,Float>)
    .def("__rdiv__", &rdiv<IFloat,Float,int,Float,Float>)
    .def("__rdiv__", &rdiv<IFloat,Float,double,Float,Float>)
    .def("__pow__", &pow<IFloat,Float,int,IFloat,int>)
    //.def("__pow__", &pow<IFloat,Float,Integer,IFloat,Integer>)
    .def("__eq__", &eq<bool,Float,double>)
    .def("__eq__", &eq<bool,Float,Float>)
    .def("__ne__", &ne<bool,Float,double>)
    .def("__ne__", &ne<bool,Float,Float>)
    .def("__lt__", &lt<bool,Float,double>)
    .def("__lt__", &lt<bool,Float,Float>)
    .def("__gt__", &gt<bool,Float,double>)
    .def("__gt__", &gt<bool,Float,Float>)
    .def("__le__", &le<bool,Float,double>)
    .def("__le__", &le<bool,Float,Float>)
    .def("__ge__", &ge<bool,Float,double>)
    .def("__ge__", &ge<bool,Float,Float>)
    //.def(self_ns::str(self))
    .def("__str__", &__str__<Float>)
    .def("__repr__", &__repr__<Float>)
  ;
  
  def("floor",&Python::floor<R,R>);
  def("ceil",&Python::ceil<R,R>);

  def("max",&Python::max<R,R,R>);
  def("min",&Python::min<R,R,R>);
  def("abs",&Python::abs<R,R>);

  def("nan",&Numeric::nan<R>);
  def("inf",&Numeric::inf<R>);
  def("infinity",&Numeric::infinity<R>);

  def("set_default_precision",&set_default_precision<Float>);
  def("default_precision",&default_precision<Float>);

}

template void export_float<Float>();
