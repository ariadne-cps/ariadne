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
#include "differentiation/multi_index.h"
#include "differentiation/taylor_variable.h"

#include "python/utilities.h"
#include "python/read_array.h"
using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

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

template<class X>
TaylorVariable<X> exp(const TaylorVariable<X>& tv) {
  return Ariadne::exp(tv);
}


template<class R>
void export_taylor_variable()
{
  typedef typename traits<R>::arithmetic_type A;
  typedef TaylorSeries<A> TS;
  typedef TaylorVariable<A> TV;


  class_<TV> taylor_variable_class(python_name<R>("TaylorVariable").c_str());
  taylor_variable_class.def("__init__", make_constructor(&make_taylor_variable<A>) );
  taylor_variable_class.def( init< uint, uint >());
  taylor_variable_class.def("__getitem__", &taylor_variable_get_item<A>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<A,double>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<A,R>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<A,A>);
  taylor_variable_class.def("__neg__", &__neg__<TV,TV>);
  taylor_variable_class.def("__add__", &__add__<TV,TV,TV>);
  taylor_variable_class.def("__add__", &__add__<TV,TV,double>);
  taylor_variable_class.def("__add__", &__add__<TV,TV,R>);
  taylor_variable_class.def("__add__", &__add__<TV,TV,A>);
  taylor_variable_class.def("__radd__", &__radd__<TV,TV,double>);
  taylor_variable_class.def("__radd__", &__radd__<TV,TV,R>);
  taylor_variable_class.def("__radd__", &__radd__<TV,TV,A>);
  taylor_variable_class.def("__sub__", &__sub__<TV,TV,TV>);
  taylor_variable_class.def("__sub__", &__sub__<TV,TV,double>);
  taylor_variable_class.def("__sub__", &__sub__<TV,TV,R>);
  taylor_variable_class.def("__sub__", &__sub__<TV,TV,A>);
  taylor_variable_class.def("__rsub__", &__rsub__<TV,TV,double>);
  taylor_variable_class.def("__rsub__", &__rsub__<TV,TV,R>);
  taylor_variable_class.def("__rsub__", &__rsub__<TV,TV,A>);
  taylor_variable_class.def("__mul__", &__mul__<TV,TV,TV>);
  taylor_variable_class.def("__mul__", &__mul__<TV,TV,double>);
  taylor_variable_class.def("__mul__", &__mul__<TV,TV,R>);
  taylor_variable_class.def("__mul__", &__mul__<TV,TV,A>);
  taylor_variable_class.def("__rmul__", &__rmul__<TV,TV,double>);
  taylor_variable_class.def("__rmul__", &__rmul__<TV,TV,R>);
  taylor_variable_class.def("__rmul__", &__rmul__<TV,TV,A>);
  taylor_variable_class.def("__div__", &__div__<TV,TV,TV>);
  taylor_variable_class.def("__div__", &__div__<TV,TV,double>);
  taylor_variable_class.def("__div__", &__div__<TV,TV,R>);
  taylor_variable_class.def("__div__", &__div__<TV,TV,A>);
  taylor_variable_class.def("__rdiv__", &__rdiv__<TV,TV,double>);
  taylor_variable_class.def("__rdiv__", &__rdiv__<TV,TV,R>);
  taylor_variable_class.def("__rdiv__", &__rdiv__<TV,TV,A>);
  taylor_variable_class.def("__pow__", &__pow__<TV,TV,int>);
  taylor_variable_class.def(self_ns::str(self));
  
  def("constant",(TV(*)(size_type, smoothness_type, const double&))&TV::constant);
  def("variable",(TV(*)(size_type, smoothness_type, const double&, size_type))&TV::variable);

  def("constant",(TV(*)(size_type, smoothness_type, const A&))&TV::constant);
  def("variable",(TV(*)(size_type, smoothness_type, const A&, size_type))&TV::variable);

  def("compose",(TV(*)(const TS&,const TV&))&compose);

  def("max",(TV(*)(const TV&,const TV&))&max<A>);
  def("min",(TV(*)(const TV&,const TV&))&min<A>);
  def("abs",(TV(*)(const TV&))&abs<A>);
  def("neg",(TV(*)(const TV&))&neg<A>);
  def("rec",(TV(*)(const TV&))&rec<A>);
  def("pow",(TV(*)(const TV&, int))&pow<A>);

  def("sqrt", (TV(*)(const TV&))&sqrt<A>);
  def("exp", (TV(*)(const TV&))&::exp<A>);
  def("log", (TV(*)(const TV&))&log<A>);
  def("sin", (TV(*)(const TV&))&sin<A>);
  def("cos", (TV(*)(const TV&))&cos<A>);
  def("tan", (TV(*)(const TV&))&tan<A>);
  def("asin", (TV(*)(const TV&))&asin<A>);
  def("acos", (TV(*)(const TV&))&acos<A>);
  def("atan", (TV(*)(const TV&))&atan<A>);

}


template<>
void export_taylor_variable<Rational>()
{
  typedef Rational Q;
  typedef TaylorSeries<Q> TS;
  typedef TaylorVariable<Q> TV;


  class_<TV> taylor_variable_class(python_name<Q>("TaylorVariable").c_str());
  taylor_variable_class.def( init< uint, uint >());
  taylor_variable_class.def("__getitem__", &taylor_variable_get_item<Q>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<Q,double>);
  taylor_variable_class.def("__setitem__",&taylor_variable_set_item<Q,Q>);
  taylor_variable_class.def("__neg__", &__neg__<TV,TV>);
  taylor_variable_class.def("__add__", &__add__<TV,TV,TV>);
  taylor_variable_class.def("__add__", &__add__<TV,TV,double>);
  taylor_variable_class.def("__add__", &__add__<TV,TV,Q>);
  taylor_variable_class.def("__radd__", &__radd__<TV,TV,double>);
  taylor_variable_class.def("__radd__", &__radd__<TV,TV,Q>);
  taylor_variable_class.def("__sub__", &__sub__<TV,TV,TV>);
  taylor_variable_class.def("__sub__", &__sub__<TV,TV,double>);
  taylor_variable_class.def("__sub__", &__sub__<TV,TV,Q>);
  taylor_variable_class.def("__rsub__", &__rsub__<TV,TV,double>);
  taylor_variable_class.def("__rsub__", &__rsub__<TV,TV,Q>);
  taylor_variable_class.def("__mul__", &__mul__<TV,TV,TV>);
  taylor_variable_class.def("__mul__", &__mul__<TV,TV,double>);
  taylor_variable_class.def("__mul__", &__mul__<TV,TV,Q>);
  taylor_variable_class.def("__rmul__", &__rmul__<TV,TV,double>);
  taylor_variable_class.def("__rmul__", &__rmul__<TV,TV,Q>);
  taylor_variable_class.def("__div__", &__div__<TV,TV,TV>);
  taylor_variable_class.def("__div__", &__div__<TV,TV,double>);
  taylor_variable_class.def("__div__", &__div__<TV,TV,Q>);
  taylor_variable_class.def("__rdiv__", &__rdiv__<TV,TV,double>);
  taylor_variable_class.def("__rdiv__", &__rdiv__<TV,TV,Q>);
  taylor_variable_class.def("__pow__", &__pow__<TV,TV,int>);
  taylor_variable_class.def(self_ns::str(self));
  
  def("constant",(TV(*)(size_type, smoothness_type, const double&))&TV::constant);
  def("variable",(TV(*)(size_type, smoothness_type, const double&, size_type))&TV::variable);

  def("constant",(TV(*)(size_type, smoothness_type, const Q&))&TV::constant);
  def("variable",(TV(*)(size_type, smoothness_type, const Q&, size_type))&TV::variable);

  def("compose",(TV(*)(const TS&,const TV&))&compose);

  def("max",(TV(*)(const TV&,const TV&))&max<Q>);
  def("min",(TV(*)(const TV&,const TV&))&min<Q>);
  def("abs",(TV(*)(const TV&))&abs<Q>);
  def("pow",(TV(*)(const TV&, int))&pow<Q>);

}




template void export_taylor_variable<Rational>();
template void export_taylor_variable<FloatPy>();
