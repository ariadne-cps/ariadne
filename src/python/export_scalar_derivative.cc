/***************************************************************************
 *            python/export_scalar_derivative.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 

#include "python/python_float.h"

#include "numeric/rational.h"
#include "numeric/interval.h"
#include "function/scalar_derivative.h"

#include <boost/python.hpp>
#include "python/python_utilities.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace Ariadne::Python;


template<class R1, class R2> inline 
void scalar_derivative_set_item(ScalarDerivative<R1>& sd, uint i, R2 x) {
  assert(i<=sd.degree()); 
  sd[i]=x;
}

template<class R> inline 
R scalar_derivative_get_item(const ScalarDerivative<R>& sd, uint i) {
  assert(i<=sd.degree()); 
  return sd[i];
}

template<class R>
void export_scalar_derivative()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef ScalarDerivative<A> SD;


  class_<SD>(python_name<R>("ScalarDerivative").c_str())
    .def( init< uint >())
    .def( init< uint, double >())
    .def( init< uint, double, double >())
    .def( init< uint, A >())
    .def( init< uint, A, A >())
    .def( init< std::string >())
    .def("__getitem__", &scalar_derivative_get_item<A>)
    .def("__setitem__",&scalar_derivative_set_item<A,double>)
    .def("__setitem__",&scalar_derivative_set_item<A,R>)
    .def("__setitem__",&scalar_derivative_set_item<A,A>)
    .def("__neg__", &Python::neg<SD,SD>)
    .def("__add__", &add<SD,SD,SD>)
    .def("__add__", &add<SD,SD,double>)
    .def("__add__", &add<SD,SD,R>)
    .def("__add__", &add<SD,SD,A>)
    .def("__radd__", &radd<SD,SD,double>)
    .def("__radd__", &radd<SD,SD,R>)
    .def("__radd__", &radd<SD,SD,A>)
    .def("__sub__", &sub<SD,SD,SD>)
    .def("__sub__", &sub<SD,SD,double>)
    .def("__sub__", &sub<SD,SD,R>)
    .def("__sub__", &sub<SD,SD,A>)
    .def("__rsub__", &rsub<SD,SD,double>)
    .def("__rsub__", &rsub<SD,SD,R>)
    .def("__rsub__", &rsub<SD,SD,A>)
    .def("__mul__", &mul<SD,SD,SD>)
    .def("__mul__", &mul<SD,SD,double>)
    .def("__mul__", &mul<SD,SD,R>)
    .def("__mul__", &mul<SD,SD,A>)
    .def("__rmul__", &rmul<SD,SD,double>)
    .def("__rmul__", &rmul<SD,SD,R>)
    .def("__rmul__", &rmul<SD,SD,A>)
    .def("__div__", &div<SD,SD,SD>)
    .def("__div__", &div<SD,SD,double>)
    .def("__div__", &div<SD,SD,R>)
    .def("__div__", &div<SD,SD,A>)
    .def("__rdiv__", &rdiv<SD,SD,double>)
    .def("__rdiv__", &rdiv<SD,SD,R>)
    .def("__rdiv__", &rdiv<SD,SD,A>)
    .def("__pow__", &pow<SD,SD,int>)
    .def(self_ns::str(self))
 ;
  
  def("compose",(SD(*)(const SD&,const SD&))&Function::compose);
  def("inverse",(SD(*)(const SD&,const A&))&Function::inverse);

  def("max",(SD(*)(const SD&,const SD&))&Function::max);
  def("min",(SD(*)(const SD&,const SD&))&Function::min);
  def("abs",(SD(*)(const SD&))&Function::abs);
  def("pow",(SD(*)(const SD&, int))&Function::pow);

  def("sqrt", (SD(*)(const SD&))&Function::sqrt);
  def("exp", (SD(*)(const SD&))&Function::exp);
  def("log", (SD(*)(const SD&))&Function::log);
  def("sin", (SD(*)(const SD&))&Function::sin);
  def("cos", (SD(*)(const SD&))&Function::cos);
  def("tan", (SD(*)(const SD&))&Function::tan);
  def("asin", (SD(*)(const SD&))&Function::asin);
  def("acos", (SD(*)(const SD&))&Function::acos);
  def("atan", (SD(*)(const SD&))&Function::atan);

}


template<>
void export_scalar_derivative<Rational>()
{
  typedef Rational Q;
  typedef ScalarDerivative<Q> SD;


  class_<SD>(python_name<Q>("ScalarDerivative").c_str())
    .def( init< std::string >())
    .def( init< uint >())
    .def( init< uint, double >())
    .def( init< uint, Q >())
    .def( init< uint, double, double >())
    .def( init< uint, Q, Q >())
    .def( init< std::string >())
    .def("__getitem__", &scalar_derivative_get_item<Q>)
    .def("__setitem__",&scalar_derivative_set_item<Q,double>)
    .def("__setitem__",&scalar_derivative_set_item<Q,Q>)
    .def("__neg__", &Python::neg<SD,SD>)
    .def("__add__", &add<SD,SD,SD>)
    .def("__add__", &add<SD,SD,double>)
    .def("__add__", &add<SD,SD,Q>)
    .def("__radd__", &radd<SD,SD,double>)
    .def("__radd__", &radd<SD,SD,Q>)
    .def("__sub__", &sub<SD,SD,SD>)
    .def("__sub__", &sub<SD,SD,double>)
    .def("__sub__", &sub<SD,SD,Q>)
    .def("__rsub__", &rsub<SD,SD,double>)
    .def("__rsub__", &rsub<SD,SD,Q>)
    .def("__mul__", &mul<SD,SD,SD>)
    .def("__mul__", &mul<SD,SD,double>)
    .def("__mul__", &mul<SD,SD,Q>)
    .def("__rmul__", &rmul<SD,SD,double>)
    .def("__rmul__", &rmul<SD,SD,Q>)
    .def("__div__", &div<SD,SD,SD>)
    .def("__div__", &div<SD,SD,double>)
    .def("__div__", &div<SD,SD,Q>)
    .def("__rdiv__", &rdiv<SD,SD,double>)
    .def("__rdiv__", &rdiv<SD,SD,Q>)
    .def("__pow__", &pow<SD,SD,int>)
    .def(self_ns::str(self))
 ;
  
  def("compose",(SD(*)(const SD&,const SD&))&Function::compose);
  def("inverse",(SD(*)(const SD&,const Q&))&Function::inverse);

  def("max",(SD(*)(const SD&,const SD&))&Function::max);
  def("min",(SD(*)(const SD&,const SD&))&Function::min);
  def("abs",(SD(*)(const SD&))&Function::abs);
  def("pow",(SD(*)(const SD&, int))&Function::pow);

}




template void export_scalar_derivative<Rational>();
template void export_scalar_derivative<Float>();
