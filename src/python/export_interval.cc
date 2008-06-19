/***************************************************************************
 *            python/export_interval.cc
 *
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

#include "python/float.h"
#include "python/name.h"
#include "python/operators.h"
#include "python/read_scalar.h"
using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R> 
std::string
__str__(const Interval<R>& ivl) {
  std::stringstream ss;
  //ss << std::fixed;
  ss << std::setprecision(20);
  ss << ivl;
  return ss.str();
}

template<class R> 
std::string
__repr__(const Interval<R>& ivl) {
  std::stringstream ss;
  //ss << std::fixed;
  ss << std::setprecision(20);
  ss << "Interval(" << ivl.lower();
  if(ivl.lower()!=ivl.upper()) {
    ss << "," << ivl.upper();
  }
  ss << ")";
  return ss.str();
}

template<class R>
Interval<R>*
make_interval(const boost::python::object& obj) 
{
  Interval<R>* ivl=new Interval<R>;
  read_scalar(*ivl,obj);
  return ivl;
}



template<class R>
void export_interval()
{
  typedef Interval<R> I;
  typedef Integer Z;

  // FIXME: pow(Interval<R>,Integer) doesn't work 

  class_< Interval<R> >(python_name<R>("Interval").c_str())
    .def("__init__", make_constructor(&make_interval<R>))
    .def(init<std::string>())
    .def(init<int,int>())
    .def(init<double,double>())
    .def(init<R,R>())
    .def(init<Rational,Rational>())
    .def(init<int>())
    .def(init<double>())
    .def(init<Rational>())
    .def(init<R>())
    .def(init< Interval<R> >())
    .def("__neg__", &__neg__<I,I>)
    .def("__add__", &__add__<I,I,int,I,R>)
    .def("__add__", &__add__<I,I,double,I,R>)
    .def("__add__", &__add__<I,I,R>)
    .def("__add__", &__add__<I,I,I>)
    .def("__radd__", &__add__<I,I,int>)
    .def("__radd__", &__add__<I,I,double>)
    .def("__radd__", &__add__<I,I,R>)
    .def("__sub__", &__sub__<I,I,int,I,R>)
    .def("__sub__", &__sub__<I,I,double,I,R>)
    .def("__sub__", &__sub__<I,I,R>)
    .def("__sub__", &__sub__<I,I,I>)
    .def("__rsub__", &__rsub__<I,I,int,I,R>)
    .def("__rsub__", &__rsub__<I,I,double,I,R>)
    .def("__rsub__", &__rsub__<I,I,R>)
    .def("__mul__", &__mul__<I,I,int,I,R>)
    .def("__mul__", &__mul__<I,I,double,I,R>)
    .def("__mul__", &__mul__<I,I,R>)
    .def("__mul__", &__mul__<I,I,I>)
    .def("__rmul__", &__mul__<I,I,int,I,R>)
    .def("__rmul__", &__mul__<I,I,double,I,R>)
    .def("__rmul__", &__mul__<I,I,R>)
    .def("__div__", &__div__<I,I,int>)
    .def("__div__", &__div__<I,I,double>)
    .def("__div__", &__div__<I,I,R>)
    .def("__div__", &__div__<I,I,I>)
    .def("__rdiv__", &__rdiv__<I,I,int,I,R>)
    .def("__rdiv__", &__rdiv__<I,I,double,I,R>)
    .def("__rdiv__", &__rdiv__<I,I,R>)
    .def("__pow__", &__pow__<I,I,int,I,int>)
    //.def("__pow__", &__pow__<I,I,Z,I,Z>)
    .def("__eq__", &__eq__<tribool,I,double>)
    .def("__eq__", &__eq__<tribool,I,R>)
    .def("__eq__", &__eq__<tribool,I,I>)
    .def("__ne__", &__ne__<tribool,I,double>)
    .def("__ne__", &__ne__<tribool,I,R>)
    .def("__ne__", &__ne__<tribool,I,I>)
    .def("__lt__", &__lt__<tribool,I,double>)
    .def("__lt__", &__lt__<tribool,I,R>)
    .def("__lt__", &__lt__<tribool,I,I>)
    .def("__gt__", &__gt__<tribool,I,double>)
    .def("__gt__", &__gt__<tribool,I,R>)
    .def("__gt__", &__gt__<tribool,I,I>)
    .def("__le__", &__le__<tribool,I,double>)
    .def("__le__", &__le__<tribool,I,R>)
    .def("__le__", &__le__<tribool,I,I>)
    .def("__ge__", &__ge__<tribool,I,double>)
    .def("__ge__", &__ge__<tribool,I,R>)
    .def("__ge__", &__ge__<tribool,I,I>)
    .def("lower", &Interval<R>::lower, return_value_policy<copy_const_reference>())
    .def("upper", &Interval<R>::upper, return_value_policy<copy_const_reference>())
    .def("midpoint", &Interval<R>::midpoint)
    .def("radius", &Interval<R>::radius)
    .def("width", &Interval<R>::width)
    .def("__str__", &__str__<FloatPy>)
    .def("__repr__", &__repr__<FloatPy>)  
  ;
 
  

  def("lower", (R(*)(const I&))&lower);
  def("upper", (R(*)(const I&))&upper);
  def("midpoint", (R(*)(const I&))&midpoint);
  def("radius", (R(*)(const I&))&radius);
  def("width", (R(*)(const I&))&radius);

  def("encloses", (bool(*)(const I&, const R&))&encloses);
  def("refines", (bool(*)(const I&, const I&))&refines);

  def("equal", (bool(*)(const I&, const I&))&equal);
  def("disjoint", (bool(*)(const I&, const I&))&disjoint);
  def("overlap", (bool(*)(const I&, const I&))&overlap);
  def("subset", (bool(*)(const I&, const I&))&subset);
  def("inside", (bool(*)(const I&, const I&))&inside);
  def("intersection", (I(*)(const I&, const I&))&intersection);
  def("hull", (I(*)(const I&, const I&))&hull);

  def("max",&Python::max<I,I,I>);
  def("min",&Python::min<I,I,I>);
  def("abs",&Python::abs<I,I>);

  def("pow",&Python::pow<I,R,int,I,int>);
  def("pow",&Python::pow<I,I,int,I,int>);

  //def("pow",&pow<I,R,Z,I,Z>);
  //def("pow",&pow<I,I,Z,I,Z>);

  def("sqrt",&sqrt<I,R,I>);
  def("sqrt",&sqrt<I,I,I>);

  def("exp",&Python::exp<I,R,I>);
  def("exp",&Python::exp<I,I,I>);
  def("log",&Python::log<I,R,I>);
  def("log",&Python::log<I,I,I>);
  def("pi",&Python::pi<I>);
  def("sin",&Python::sin<I,R,I>);
  def("sin",&Python::sin<I,I,I>);
  def("cos",&Python::cos<I,R,I>);
  def("cos",&Python::cos<I,I,I>);
  def("tan",&Python::tan<I,R,I>);
  def("tan",&Python::tan<I,I,I>);
  def("asin",&Python::asin<I,R,I>);
  def("asin",&Python::asin<I,I,I>);
  def("acos",&Python::acos<I,R,I>);
  def("acos",&Python::acos<I,I,I>);
  def("atan",&Python::atan<I,R,I>);
  def("atan",&Python::atan<I,I,I>);

}

template void export_interval<FloatPy>();
