/***************************************************************************
 *            python/export_taylor_derivative.cc
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
#include "function/multi_index.h"
#include "function/taylor_derivative.h"

#include <boost/python.hpp>
#include "python/python_utilities.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace Ariadne::Python;


template<class R1, class R2> inline 
void taylor_derivative_set_item(TaylorDerivative<R1>& td, const MultiIndex& i, R2 x) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  td[i]=x;
}

template<class R> inline 
R taylor_derivative_get_item(const TaylorDerivative<R>& td, const MultiIndex& i) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  return td[i];
}

template<class R>
void export_taylor_derivative()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef ScalarDerivative<A> SD;
  typedef TaylorDerivative<A> TD;


  class_<TD>(python_name<R>("TaylorDerivative").c_str())
    .def( init< uint, uint >())
    .def( init< uint, uint, double >())
    .def( init< uint, uint, uint, double >())
    .def( init< uint, uint, A >())
    .def( init< uint, uint, uint, A >())
    .def("__getitem__", &taylor_derivative_get_item<A>)
    .def("__setitem__",&taylor_derivative_set_item<A,double>)
    .def("__setitem__",&taylor_derivative_set_item<A,R>)
    .def("__setitem__",&taylor_derivative_set_item<A,A>)
    .def("__neg__", &Python::neg<TD,TD>)
    .def("__add__", &add<TD,TD,TD>)
    .def("__add__", &add<TD,TD,double>)
    .def("__add__", &add<TD,TD,R>)
    .def("__add__", &add<TD,TD,A>)
    .def("__radd__", &radd<TD,TD,double>)
    .def("__radd__", &radd<TD,TD,R>)
    .def("__radd__", &radd<TD,TD,A>)
    .def("__sub__", &sub<TD,TD,TD>)
    .def("__sub__", &sub<TD,TD,double>)
    .def("__sub__", &sub<TD,TD,R>)
    .def("__sub__", &sub<TD,TD,A>)
    .def("__rsub__", &rsub<TD,TD,double>)
    .def("__rsub__", &rsub<TD,TD,R>)
    .def("__rsub__", &rsub<TD,TD,A>)
    .def("__mul__", &mul<TD,TD,TD>)
    .def("__mul__", &mul<TD,TD,double>)
    .def("__mul__", &mul<TD,TD,R>)
    .def("__mul__", &mul<TD,TD,A>)
    .def("__rmul__", &rmul<TD,TD,double>)
    .def("__rmul__", &rmul<TD,TD,R>)
    .def("__rmul__", &rmul<TD,TD,A>)
    .def("__div__", &div<TD,TD,TD>)
    .def("__div__", &div<TD,TD,double>)
    .def("__div__", &div<TD,TD,R>)
    .def("__div__", &div<TD,TD,A>)
    .def("__rdiv__", &rdiv<TD,TD,double>)
    .def("__rdiv__", &rdiv<TD,TD,R>)
    .def("__rdiv__", &rdiv<TD,TD,A>)
    .def("__pow__", &pow<TD,TD,int>)
    .def(self_ns::str(self))
 ;
  
  def("compose",(TD(*)(const SD&,const TD&))&Function::compose);

  def("max",(TD(*)(const TD&,const TD&))&Function::max);
  def("min",(TD(*)(const TD&,const TD&))&Function::min);
  def("abs",(TD(*)(const TD&))&Function::abs);
  def("pow",(TD(*)(const TD&, int))&Function::pow);

  def("sqrt", (TD(*)(const TD&))&Function::sqrt);
  def("exp", (TD(*)(const TD&))&Function::exp);
  def("log", (TD(*)(const TD&))&Function::log);
  def("sin", (TD(*)(const TD&))&Function::sin);
  def("cos", (TD(*)(const TD&))&Function::cos);
  def("tan", (TD(*)(const TD&))&Function::tan);
  def("asin", (TD(*)(const TD&))&Function::asin);
  def("acos", (TD(*)(const TD&))&Function::acos);
  def("atan", (TD(*)(const TD&))&Function::atan);

}


template<>
void export_taylor_derivative<Rational>()
{
  typedef Rational Q;
  typedef ScalarDerivative<Q> SD;
  typedef TaylorDerivative<Q> TD;


  class_<TD>(python_name<Q>("TaylorDerivative").c_str())
    .def( init< uint, uint >())
    .def( init< uint, uint, double >())
    .def( init< uint, uint, Q >())
    .def( init< uint, uint, uint, double >())
    .def( init< uint, uint, uint, Q >())
    .def("__getitem__", &taylor_derivative_get_item<Q>)
    .def("__setitem__",&taylor_derivative_set_item<Q,double>)
    .def("__setitem__",&taylor_derivative_set_item<Q,Q>)
    .def("__neg__", &Python::neg<TD,TD>)
    .def("__add__", &add<TD,TD,TD>)
    .def("__add__", &add<TD,TD,double>)
    .def("__add__", &add<TD,TD,Q>)
    .def("__radd__", &radd<TD,TD,double>)
    .def("__radd__", &radd<TD,TD,Q>)
    .def("__sub__", &sub<TD,TD,TD>)
    .def("__sub__", &sub<TD,TD,double>)
    .def("__sub__", &sub<TD,TD,Q>)
    .def("__rsub__", &rsub<TD,TD,double>)
    .def("__rsub__", &rsub<TD,TD,Q>)
    .def("__mul__", &mul<TD,TD,TD>)
    .def("__mul__", &mul<TD,TD,double>)
    .def("__mul__", &mul<TD,TD,Q>)
    .def("__rmul__", &rmul<TD,TD,double>)
    .def("__rmul__", &rmul<TD,TD,Q>)
    .def("__div__", &div<TD,TD,TD>)
    .def("__div__", &div<TD,TD,double>)
    .def("__div__", &div<TD,TD,Q>)
    .def("__rdiv__", &rdiv<TD,TD,double>)
    .def("__rdiv__", &rdiv<TD,TD,Q>)
    .def("__pow__", &pow<TD,TD,int>)
    .def(self_ns::str(self))
 ;
  
  def("compose",(TD(*)(const SD&,const TD&))&Function::compose);

  def("max",(TD(*)(const TD&,const TD&))&Function::max);
  def("min",(TD(*)(const TD&,const TD&))&Function::min);
  def("abs",(TD(*)(const TD&))&Function::abs);
  def("pow",(TD(*)(const TD&, int))&Function::pow);

}




template void export_taylor_derivative<Rational>();
template void export_taylor_derivative<Float>();
