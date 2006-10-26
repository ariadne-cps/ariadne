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

using namespace boost::python;
using namespace Ariadne;
using Numeric::Interval;

template<class R> inline Interval<R> neginvl(const Interval<R>& i) { return -i; }
template<class R> inline Interval<R> addinvl(const Interval<R>& i, const Interval<R>& j) { return i+j; }
template<class R> inline Interval<R> subinvl(const Interval<R>& i, const Interval<R>& j) { return i-j; }
template<class R> inline Interval<R> mulinvl(const Interval<R>& i, const Interval<R>& j) { return i*j; }
template<class R> inline Interval<R> divinvl(const Interval<R>& i, const Interval<R>& j) { return i/j; }

template<class R1, class R2>
inline Interval<R1> addreal(const Interval<R1>& i, const R2& r) { return i+(R1)r; }
template<class R1, class R2>
inline Interval<R1> subreal(const Interval<R1>& i, const R2& r) { return i-(R1)r; }
template<class R1, class R2>
inline Interval<R1> mulreal(const Interval<R1>& i, const R2& r) { return i*(R1)r; }
template<class R1, class R2>
inline Interval<R1> divreal(const Interval<R1>& i, const R2& r) { return i/(R1)r; }

template<class R>
void export_interval()
{
  class_< Interval<R> >(python_name<R>("Interval").c_str())
    .def(init<int,int>())
    .def(init<double,double>())
    .def(init<R,R>())
    .def(init<int>())
    .def(init<double>())
    .def(init<R>())
    .def("__neg__", &neginvl<R>)
    .def("__add__", &addinvl<R>) 
    .def("__sub__", &subinvl<R>) 
    .def("__mul__", &mulinvl<R>) 
    .def("__div__", &divinvl<R>)  // __div__
    .def("__add__", &addreal<R,int>) 
    .def("__sub__", &subreal<R,int>) 
    .def("__mul__", &mulreal<R,int>) 
    .def("__div__", &divreal<R,int>)
    .def("__add__", &addreal<R,float>) 
    .def("__sub__", &subreal<R,float>) 
    .def("__mul__", &mulreal<R,float>) 
    .def("__div__", &divreal<R,float>)
    .def("__add__", &addreal<R,R>) 
    .def("__sub__", &subreal<R,R>) 
    .def("__mul__", &mulreal<R,R>) 
    .def("__div__", &divreal<R,R>)
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

template void export_interval<Float64>();
template void export_interval<MPFloat>();
