/***************************************************************************
 *            python/export_taylor_variable.cc
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
 

#include "python/float.h"

#include "numeric/rational.h"
#include "numeric/interval.h"
#include "function/multi_index.h"
#include "function/taylor_variable.h"

#include <boost/python.hpp>
#include "python/utilities.h"
#include "python/read_array.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace Ariadne::Python;

template<class X>
TaylorVariable<X>*
make_taylor_variable(const uint& as, const uint& d, const boost::python::object& obj) 
{
  TaylorVariable<X>* result=new TaylorVariable<X>(as,d);
  read_array(result->data(),obj);
  assert(result->data().size()==compute_polynomial_data_size(as,d));
  return result;
}


template<class R1, class R2> inline 
void taylor_variable_set_item(TaylorVariable<R1>& td, const MultiIndex& i, R2 x) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  td[i]=x;
}

template<class R> inline 
R taylor_variable_get_item(const TaylorVariable<R>& td, const MultiIndex& i) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  return td[i];
}

template<class R>
void export_taylor_variable()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef TaylorSeries<A> TS;
  typedef TaylorVariable<A> TD;


  class_<TD> taylor_variable_class(python_name<R>("TaylorVariable").c_str());
  taylor_variable_class.def("__init__", make_constructor(&make_taylor_variable<A>) );
  taylor_variable_class.def( init< uint, uint >());
  taylor_variable_class.def("__getitem__", &taylor_variable_get_item<A>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<A,double>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<A,R>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<A,A>);
  taylor_variable_class.def("__neg__", &Python::neg<TD,TD>);
  taylor_variable_class.def("__add__", &Python::add<TD,TD,TD>);
  taylor_variable_class.def("__add__", &Python::add<TD,TD,double>);
  taylor_variable_class.def("__add__", &Python::add<TD,TD,R>);
  taylor_variable_class.def("__add__", &Python::add<TD,TD,A>);
  taylor_variable_class.def("__radd__", &Python::radd<TD,TD,double>);
  taylor_variable_class.def("__radd__", &Python::radd<TD,TD,R>);
  taylor_variable_class.def("__radd__", &Python::radd<TD,TD,A>);
  taylor_variable_class.def("__sub__", &Python::sub<TD,TD,TD>);
  taylor_variable_class.def("__sub__", &Python::sub<TD,TD,double>);
  taylor_variable_class.def("__sub__", &Python::sub<TD,TD,R>);
  taylor_variable_class.def("__sub__", &Python::sub<TD,TD,A>);
  taylor_variable_class.def("__rsub__", &Python::rsub<TD,TD,double>);
  taylor_variable_class.def("__rsub__", &Python::rsub<TD,TD,R>);
  taylor_variable_class.def("__rsub__", &Python::rsub<TD,TD,A>);
  taylor_variable_class.def("__mul__", &Python::mul<TD,TD,TD>);
  taylor_variable_class.def("__mul__", &Python::mul<TD,TD,double>);
  taylor_variable_class.def("__mul__", &Python::mul<TD,TD,R>);
  taylor_variable_class.def("__mul__", &Python::mul<TD,TD,A>);
  taylor_variable_class.def("__rmul__", &Python::rmul<TD,TD,double>);
  taylor_variable_class.def("__rmul__", &Python::rmul<TD,TD,R>);
  taylor_variable_class.def("__rmul__", &Python::rmul<TD,TD,A>);
  taylor_variable_class.def("__div__", &Python::div<TD,TD,TD>);
  taylor_variable_class.def("__div__", &Python::div<TD,TD,double>);
  taylor_variable_class.def("__div__", &Python::div<TD,TD,R>);
  taylor_variable_class.def("__div__", &Python::div<TD,TD,A>);
  taylor_variable_class.def("__rdiv__", &Python::rdiv<TD,TD,double>);
  taylor_variable_class.def("__rdiv__", &Python::rdiv<TD,TD,R>);
  taylor_variable_class.def("__rdiv__", &Python::rdiv<TD,TD,A>);
  taylor_variable_class.def("__pow__", &Python::pow<TD,TD,int>);
  taylor_variable_class.def(self_ns::str(self));
  
  def("compose",(TD(*)(const TS&,const TD&))&Function::compose);

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
void export_taylor_variable<Rational>()
{
  typedef Rational Q;
  typedef TaylorSeries<Q> TS;
  typedef TaylorVariable<Q> TD;


  class_<TD> taylor_variable_class(python_name<Q>("TaylorVariable").c_str());
  taylor_variable_class.def( init< uint, uint >());
  taylor_variable_class.def("__getitem__", &taylor_variable_get_item<Q>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<Q,double>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<Q,Q>);
  taylor_variable_class.def("__neg__", &Python::neg<TD,TD>);
  taylor_variable_class.def("__add__", &Python::add<TD,TD,TD>);
  taylor_variable_class.def("__add__", &Python::add<TD,TD,double>);
  taylor_variable_class.def("__add__", &Python::add<TD,TD,Q>);
  taylor_variable_class.def("__radd__", &Python::radd<TD,TD,double>);
  taylor_variable_class.def("__radd__", &Python::radd<TD,TD,Q>);
  taylor_variable_class.def("__sub__", &Python::sub<TD,TD,TD>);
  taylor_variable_class.def("__sub__", &Python::sub<TD,TD,double>);
  taylor_variable_class.def("__sub__", &Python::sub<TD,TD,Q>);
  taylor_variable_class.def("__rsub__", &Python::rsub<TD,TD,double>);
  taylor_variable_class.def("__rsub__", &Python::rsub<TD,TD,Q>);
  taylor_variable_class.def("__mul__", &Python::mul<TD,TD,TD>);
  taylor_variable_class.def("__mul__", &Python::mul<TD,TD,double>);
  taylor_variable_class.def("__mul__", &Python::mul<TD,TD,Q>);
  taylor_variable_class.def("__rmul__", &Python::rmul<TD,TD,double>);
  taylor_variable_class.def("__rmul__", &Python::rmul<TD,TD,Q>);
  taylor_variable_class.def("__div__", &Python::div<TD,TD,TD>);
  taylor_variable_class.def("__div__", &Python::div<TD,TD,double>);
  taylor_variable_class.def("__div__", &Python::div<TD,TD,Q>);
  taylor_variable_class.def("__rdiv__", &Python::rdiv<TD,TD,double>);
  taylor_variable_class.def("__rdiv__", &Python::rdiv<TD,TD,Q>);
  taylor_variable_class.def("__pow__", &Python::pow<TD,TD,int>);
  taylor_variable_class.def(self_ns::str(self));
  
  def("compose",(TD(*)(const TS&,const TD&))&Function::compose);

  def("max",(TD(*)(const TD&,const TD&))&Function::max);
  def("min",(TD(*)(const TD&,const TD&))&Function::min);
  def("abs",(TD(*)(const TD&))&Function::abs);
  def("pow",(TD(*)(const TD&, int))&Function::pow);

}




template void export_taylor_variable<Rational>();
template void export_taylor_variable<FloatPy>();
