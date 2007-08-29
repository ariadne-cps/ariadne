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
using namespace Ariadne::Python;


template<class R>
void export_interval()
{
  typedef Interval<R> I;
  typedef Integer Z;

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
    .def("__pow__", &pow<I,I,int,I,int>)
    .def("__pow__", &pow<I,I,Z,I,Z>)
    .def("__eq__", &eq<tribool,I,double>)
    .def("__eq__", &eq<tribool,I,R>)
    .def("__eq__", &eq<tribool,I,I>)
    .def("__ne__", &ne<tribool,I,double>)
    .def("__ne__", &ne<tribool,I,R>)
    .def("__ne__", &ne<tribool,I,I>)
    .def("__lt__", &lt<tribool,I,double>)
    .def("__lt__", &lt<tribool,I,R>)
    .def("__lt__", &lt<tribool,I,I>)
    .def("__gt__", &gt<tribool,I,double>)
    .def("__gt__", &gt<tribool,I,R>)
    .def("__gt__", &gt<tribool,I,I>)
    .def("__le__", &le<tribool,I,double>)
    .def("__le__", &le<tribool,I,R>)
    .def("__le__", &le<tribool,I,I>)
    .def("__ge__", &ge<tribool,I,double>)
    .def("__ge__", &ge<tribool,I,R>)
    .def("__ge__", &ge<tribool,I,I>)
    .def("lower", &Interval<R>::lower, return_value_policy<copy_const_reference>())
    .def("upper", &Interval<R>::upper, return_value_policy<copy_const_reference>())
    .def("midpoint", &Interval<R>::midpoint)
    .def("radius", &Interval<R>::radius)
    .def("width", &Interval<R>::width)
    .def(self_ns::str(self))    // __self_ns::str__
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

  def("pow",&pow<I,R,int,I,int>);
  def("pow",&pow<I,I,int,I,int>);

  def("pow",&pow<I,R,Z,I,Z>);
  def("pow",&pow<I,I,Z,I,Z>);

  def("sqrt",&sqrt<I,R,I>);
  def("sqrt",&sqrt<I,I,I>);

  def("exp",&exp<I,R,I>);
  def("exp",&exp<I,I,I>);
  def("log",&log<I,R,I>);
  def("log",&log<I,I,I>);
  def("pi",&Numeric::pi<R>);
  def("sin",&sin<I,R,I>);
  def("sin",&sin<I,I,I>);
  def("cos",&cos<I,R,I>);
  def("cos",&cos<I,I,I>);
  def("tan",&tan<I,R,I>);
  def("tan",&tan<I,I,I>);
  def("asin",&asin<I,R,I>);
  def("asin",&asin<I,I,I>);
  def("acos",&acos<I,R,I>);
  def("acos",&acos<I,I,I>);
  def("atan",&atan<I,R,I>);
  def("atan",&atan<I,I,I>);

}

template void export_interval<Float>();
