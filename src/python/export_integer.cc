/***************************************************************************
 *            python/export_integer.cc
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
#include "numeric/integer.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

std::string
__str__(const Integer& z)
{
  std::stringstream ss;
  ss << z;
  return ss.str();
}

std::string
__repr__(const Integer& z)
{
  std::stringstream ss;
  ss << "Integer(" << z << ")";
  return ss.str();
}


void export_integer() {
  class_<Integer>("Integer")
    .def(init<int>())
    .def(init<std::string>())
    .def(init<Integer>())
    .def("__neg__", &__neg__<Integer,Integer>)
    .def("__add__", &__add__<Integer,Integer,Integer>)
    .def("__add__", &__add__<Integer,Integer,int>)
    .def("__radd__", &__add__<Integer,Integer,int>)
    .def("__sub__", &__sub__<Integer,Integer,Integer>)
    .def("__sub__", &__sub__<Integer,Integer,int>)
    .def("__rsub__", &__rsub__<Integer,Integer,int>)
    .def("__mul__", &__mul__<Integer,Integer,Integer>)
    .def("__mul__", &__mul__<Integer,Integer,int>)
    .def("__rmul__", &__mul__<Integer,Integer,int>)
    .def("__div__", &__div__<Rational,Integer,Integer,Rational,Rational>)
    .def("__div__", &__div__<Rational,Integer,int,Rational,Rational>)
    .def("__rdiv__", &__rdiv__<Rational,Integer,int,Rational,Rational>)
    .def("__eq__", &__eq__<bool,Integer,int>)
    .def("__eq__", &__eq__<bool,Integer,Integer>)
    .def("__ne__", &__ne__<bool,Integer,int>)
    .def("__ne__", &__ne__<bool,Integer,Integer>)
    .def("__lt__", &__lt__<bool,Integer,int>)
    .def("__lt__", &__lt__<bool,Integer,Integer>)
    .def("__gt__", &__gt__<bool,Integer,int>)
    .def("__gt__", &__gt__<bool,Integer,Integer>)
    .def("__le__", &__le__<bool,Integer,int>)
    .def("__le__", &__le__<bool,Integer,Integer>)
    .def("__ge__", &__ge__<bool,Integer,int>)
    .def("__ge__", &__ge__<bool,Integer,Integer>)
    .def("__str__", &__str__)
    .def("__repr__", &__repr__)
  ;

  def("max",&Python::max<Integer,Integer,Integer>);
  def("min",&Python::min<Integer,Integer,Integer>);
  def("abs",&Python::abs<Integer,Integer>);

  def("pow",&Python::pow<Integer,Integer,uint>);
  
}
