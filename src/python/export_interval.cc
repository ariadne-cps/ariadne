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

#include "python/typedefs.h"

using namespace boost::python;
using namespace Ariadne;

typedef RInterval (*IntervalFunc) (const RInterval&);
typedef RInterval (*IntervalBinFunc) (const RInterval&, const RInterval&);

/* FIXME: Do these divisions */
inline RInterval div_ii(const RInterval& r, const RInterval& i) { assert(false); }
inline RInterval div_ri(const Real& r, const RInterval& i) { assert(false); }

inline RInterval neginvl(const RInterval& i) { return -i; }
inline RInterval addinvl(const RInterval& i, const RInterval& j) { return i+j; }
inline RInterval subinvl(const RInterval& i, const RInterval& j) { return i-j; }
inline RInterval mulinvl(const RInterval& i, const RInterval& j) { return i*j; }

template <typename R>
inline RInterval addreal(const RInterval& i, const R& r) { return i+(Real)r; }
template <typename R>
inline RInterval subreal(const RInterval& i, const R& r) { return i-(Real)r; }
template <typename R>
inline RInterval mulreal(const RInterval& i, const R& r) { return i*(Real)r; }
template <typename R>
inline RInterval div_ir(const RInterval& i, const R& r) { return i/(Real)r; }

void export_interval()
{
  class_< RInterval >("Interval")
    .def(init<int,int>())
    .def(init<double,double>())
    .def(init<Real,Real>())
//    .def(init<double>())
    .def(init<Real>())
    .def("__neg__", &neginvl)
    .def("__add__", &addinvl) 
    .def("__sub__", &subinvl) 
    .def("__mul__", &mulinvl) 
    .def("__div__",&div_ii)  // __div__
    .def("__add__", &addreal<int>) 
    .def("__sub__", &subreal<int>) 
    .def("__mul__", &mulreal<int>) 
    .def("__div__",&div_ir<int>)
    .def("__add__", &addreal<float>) 
    .def("__sub__", &subreal<float>) 
    .def("__mul__", &mulreal<float>) 
    .def("__div__",&div_ir<float>)
    .def("__add__", &addreal<Real>) 
    .def("__sub__", &subreal<Real>) 
    .def("__mul__", &mulreal<Real>) 
    .def("__div__",&div_ir<Real>)
    .def("__div__",&div_ri)  // __div__
    .def("lower", &RInterval::lower, return_value_policy<copy_const_reference>())
    .def("upper", &RInterval::upper, return_value_policy<copy_const_reference>())
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
