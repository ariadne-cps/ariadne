/***************************************************************************
 *            python/export_rational.cc
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


#include "python/operators.h"
#include "python/float.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

std::string
__str__(const Rational& q)
{
  std::stringstream ss;
  ss << q;
  return ss.str();
}

std::string
__repr__(const Rational& q)
{
  std::stringstream ss;
  ss << "Rational(" << q.numerator();
  if(q!=q.numerator()) {
    ss << "," << q.denominator();
  }
  ss << ")";
  return ss.str();
}


void export_rational() {
  class_<Rational>("Rational")
    .def(init<std::string>())
    .def(init<int,int>())
    .def(init<Integer,Integer>())
    .def(init<int>())
    .def(init<Integer>())
    .def(init<double>())
    .def(init<FloatPy>())
    .def(init<Rational>())
    .def("__neg__", &__neg__<Rational,Rational>)
    .def("__add__", &__add__<Rational,Rational,int>)
    .def("__add__", &__add__<Rational,Rational,Integer>)
    .def("__add__", &__add__<Rational,Rational,double>)
    .def("__add__", &__add__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__add__", &__add__<Rational,Rational,Rational>)
    .def("__radd__", &__add__<Rational,Rational,int>)
    .def("__radd__", &__add__<Rational,Rational,Integer>)
    .def("__radd__", &__add__<Rational,Rational,double>)
    .def("__radd__", &__add__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__sub__", &__sub__<Rational,Rational,int>)
    .def("__sub__", &__sub__<Rational,Rational,Integer>)
    .def("__sub__", &__sub__<Rational,Rational,double>)
    .def("__sub__", &__sub__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__sub__", &__sub__<Rational,Rational,Rational>)
    .def("__rsub__", &__rsub__<Rational,Rational,int>)
    .def("__rsub__", &__rsub__<Rational,Rational,Integer>)
    .def("__rsub__", &__rsub__<Rational,Rational,double>)
    .def("__rsub__", &__rsub__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__mul__", &__mul__<Rational,Rational,int>)
    .def("__mul__", &__mul__<Rational,Rational,Integer>)
    .def("__mul__", &__mul__<Rational,Rational,double>)
    .def("__mul__", &__mul__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__mul__", &__mul__<Rational,Rational,Rational>)
    .def("__rmul__", &__mul__<Rational,Rational,int>)
    .def("__rmul__", &__mul__<Rational,Rational,Integer>)
    .def("__rmul__", &__mul__<Rational,Rational,double>)
    .def("__rmul__", &__mul__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__div__", &__div__<Rational,Rational,int>)
    .def("__div__", &__div__<Rational,Rational,Integer>)
    .def("__div__", &__div__<Rational,Rational,double>)
    .def("__div__", &__div__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__div__", &__div__<Rational,Rational,Rational>)
    .def("__rdiv__", &__rdiv__<Rational,Rational,int>)
    .def("__rdiv__", &__rdiv__<Rational,Rational,Integer>)
    .def("__rdiv__", &__rdiv__<Rational,Rational,double>)
    .def("__rdiv__", &__rdiv__<Rational,Rational,FloatPy,Rational,Rational>)
    .def("__eq__", &__eq__<bool,Rational,Rational>)
    .def("__eq__", &__eq__<bool,Rational,double>)
    .def("__eq__", &__eq__<bool,Rational,Integer>)
    .def("__eq__", &__eq__<bool,Rational,double>)
    .def("__ne__", &__ne__<bool,Rational,Rational>)
    .def("__ne__", &__ne__<bool,Rational,double>)
    .def("__ne__", &__ne__<bool,Rational,Integer>)
    .def("__ne__", &__ne__<bool,Rational,double>)
    .def("__lt__", &__lt__<bool,Rational,Rational>)
    .def("__lt__", &__lt__<bool,Rational,double>)
    .def("__lt__", &__lt__<bool,Rational,Integer>)
    .def("__lt__", &__lt__<bool,Rational,double>)
    .def("__gt__", &__gt__<bool,Rational,double>)
    .def("__gt__", &__gt__<bool,Rational,Integer>)
    .def("__gt__", &__gt__<bool,Rational,double>)
    .def("__le__", &__le__<bool,Rational,Rational>)
    .def("__le__", &__le__<bool,Rational,double>)
    .def("__le__", &__le__<bool,Rational,Integer>)
    .def("__le__", &__le__<bool,Rational,double>)
    .def("__ge__", &__ge__<bool,Rational,double>)
    .def("__ge__", &__ge__<bool,Rational,Integer>)
    .def("__ge__", &__ge__<bool,Rational,double>)
    .def("numerator", &Rational::numerator)
    .def("denominator", &Rational::denominator)
    .def("__str__", &__str__)
    .def("__repr__", &__repr__)
  ;
 
  def("max",&Python::max<Rational,Rational,Rational>);
  def("min",&Python::min<Rational,Rational,Rational>);
  def("abs",&Python::abs<Rational,Rational>);
  def("floor",&Python::floor<Integer,Rational>);
  def("ceil",&Python::ceil<Integer,Rational>);
    
}
