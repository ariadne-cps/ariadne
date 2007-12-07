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

#include "python/utilities.h"

#include "python/float.h"
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
  if(x==Numeric::floor<R>(x)) {
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
  ss << "FloatPy(";
  if(x==Numeric::floor<R>(x)) {
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
  typedef R FloatPy;
  typedef Interval<R> IFloat;
  
  def("set_output_precision",(void(*)(uint))&set_output_precision);

  class_<FloatPy>("Float")
    //.def("__init__", make_constructor(&make_float<R>) )
    .def(init<int>())
    .def(init<double>())
    .def(init<std::string>())
    .def(init<FloatPy>())
    .def("__neg__", &neg<FloatPy,FloatPy>)
    .def("__add__", &add<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__add__", &add<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__add__", &add<IFloat,FloatPy,FloatPy>)
    .def("__radd__", &add<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__radd__", &add<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__sub__", &sub<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__sub__", &sub<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__sub__", &sub<IFloat,FloatPy,FloatPy>)
    .def("__rsub__", &rsub<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__rsub__", &rsub<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__mul__", &mul<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__mul__", &mul<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__mul__", &mul<IFloat,FloatPy,FloatPy>)
    .def("__rmul__", &mul<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__rmul__", &mul<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__div__", &div<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__div__", &div<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__div__", &div<IFloat,FloatPy,FloatPy>)
    .def("__rdiv__", &rdiv<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__rdiv__", &rdiv<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__pow__", &pow<IFloat,FloatPy,int,IFloat,int>)
    //.def("__pow__", &pow<IFloat,FloatPy,Integer,IFloat,Integer>)
    .def("__eq__", &eq<bool,FloatPy,double>)
    .def("__eq__", &eq<bool,FloatPy,FloatPy>)
    .def("__ne__", &ne<bool,FloatPy,double>)
    .def("__ne__", &ne<bool,FloatPy,FloatPy>)
    .def("__lt__", &lt<bool,FloatPy,double>)
    .def("__lt__", &lt<bool,FloatPy,FloatPy>)
    .def("__gt__", &gt<bool,FloatPy,double>)
    .def("__gt__", &gt<bool,FloatPy,FloatPy>)
    .def("__le__", &le<bool,FloatPy,double>)
    .def("__le__", &le<bool,FloatPy,FloatPy>)
    .def("__ge__", &ge<bool,FloatPy,double>)
    .def("__ge__", &ge<bool,FloatPy,FloatPy>)
    //.def(self_ns::str(self))
    .def("__str__", &__str__<FloatPy>)
    .def("__repr__", &__repr__<FloatPy>)
  ;
  
  def("floor",&Python::floor<R,R>);
  def("ceil",&Python::ceil<R,R>);

  def("max",&Python::max<R,R,R>);
  def("min",&Python::min<R,R,R>);
  def("abs",&Python::abs<R,R>);

  def("nan",&Numeric::nan<R>);
  def("inf",&Numeric::inf<R>);

  def("set_default_precision",&set_default_precision<FloatPy>);
  def("default_precision",&default_precision<FloatPy>);

}

template void export_float<FloatPy>();
