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

#include "python/operators.h"

#include "python/float.h"
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/interval.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

void set_output_precision(uint);
template<class R> void set_default_precision(uint);
template<class R> uint default_precision();

void set_output_precision(uint p) {
  std::cout << std::setprecision(p);
}


#if defined PYTHON_FLOAT64

#warning "Compiling with FloatPy == Float64"

template<> void set_default_precision<FloatPy>(uint p) { 
  throw std::runtime_error("Cannot set precision of 64-bit float");
}

template<> uint default_precision<FloatPy>() { 
  return 64;
}

#elif defined PYTHON_FLOATMP

#warning "Compiling with FloatPy == FloatMP"

template<> void set_default_precision<FloatPy>(uint p) { 
  FloatPy::set_default_precision(p);
}

template<> uint default_precision<FloatPy>() { 
  return FloatPy::default_precision();
}

#endif

template<class R> 
std::string
__str__(const R& x) {
  std::stringstream ss;
  if(x==Python::floor<R>(x)) {
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
  if(x==Python::floor<R>(x)) {
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

  class_<FloatPy>("Float",init<>())
    //.def("__init__", make_constructor(&make_float<R>) )
    .def(init<int>())
    .def(init<double>())
    .def(init<std::string>())
    .def(init<FloatPy>())
    .def("get_d", &FloatPy::get_d)
    .def("__neg__", &__neg__<FloatPy,FloatPy>)
    .def("__add__", &__add__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__add__", &__add__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__add__", &__add__<IFloat,FloatPy,FloatPy>)
    .def("__radd__", &__add__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__radd__", &__add__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__sub__", &__sub__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__sub__", &__sub__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__sub__", &__sub__<IFloat,FloatPy,FloatPy>)
    .def("__rsub__", &__rsub__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__rsub__", &__rsub__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__mul__", &__mul__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__mul__", &__mul__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__mul__", &__mul__<IFloat,FloatPy,FloatPy>)
    .def("__rmul__", &__mul__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__rmul__", &__mul__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__div__", &__div__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__div__", &__div__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__div__", &__div__<IFloat,FloatPy,FloatPy>)
    .def("__rdiv__", &__rdiv__<IFloat,FloatPy,int,FloatPy,FloatPy>)
    .def("__rdiv__", &__rdiv__<IFloat,FloatPy,double,FloatPy,FloatPy>)
    .def("__pow__", &__pow__<IFloat,FloatPy,int,IFloat,int>)
    //.def("__pow__", &__pow__<IFloat,FloatPy,Integer,IFloat,Integer>)
    .def("__eq__", &__eq__<bool,FloatPy,double>)
    .def("__eq__", &__eq__<bool,FloatPy,FloatPy>)
    .def("__ne__", &__ne__<bool,FloatPy,double>)
    .def("__ne__", &__ne__<bool,FloatPy,FloatPy>)
    .def("__lt__", &__lt__<bool,FloatPy,double>)
    .def("__lt__", &__lt__<bool,FloatPy,FloatPy>)
    .def("__gt__", &__gt__<bool,FloatPy,double>)
    .def("__gt__", &__gt__<bool,FloatPy,FloatPy>)
    .def("__le__", &__le__<bool,FloatPy,double>)
    .def("__le__", &__le__<bool,FloatPy,FloatPy>)
    .def("__ge__", &__ge__<bool,FloatPy,double>)
    .def("__ge__", &__ge__<bool,FloatPy,FloatPy>)
    //.def(self_ns::str(self))
    .def("__str__", &__str__<FloatPy>)
    .def("__repr__", &__repr__<FloatPy>)
  ;
  
  def("floor",&Python::floor<R,R>);
  def("ceil",&Python::ceil<R,R>);

  def("max",&Python::max<R,R,R>);
  def("min",&Python::min<R,R,R>);
  def("abs",&Python::abs<R,R>);

  def("nan",&nan<R>);
  def("inf",&inf<R>);

  def("set_default_precision",&set_default_precision<FloatPy>);
  def("default_precision",&default_precision<FloatPy>);

}

template void export_float<FloatPy>();
