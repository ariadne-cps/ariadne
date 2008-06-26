/***************************************************************************
 *            python/export_taylor_series.cc
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
#include "differentiation/taylor_series.h"
#include "differentiation/taylor_variable.h"
#include "python/read_array.h"
using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include "python/utilities.h"
using namespace boost::python;


template<class X>
TaylorSeries<X>*
make_taylor_series(const boost::python::object& obj) 
{
  TaylorSeries<X>* ts=new TaylorSeries<X>;
  read_array(ts->data(),obj);
  return ts;
}

template<class R1, class R2> inline 
void set_item(TaylorSeries<R1>& sd, uint i, R2 x) {
  assert(i<=sd.degree()); 
  sd[i]=x;
}

template<class R> inline 
R get_item(const TaylorSeries<R>& sd, uint i) {
  assert(i<=sd.degree()); 
  return sd[i];
}

template<class X> inline
TaylorSeries<X> exp(const TaylorSeries<X>& ts) {
  return Ariadne::exp(ts);
}


template<class R>
void export_taylor_series()
{
  typedef typename traits<R>::arithmetic_type A;
  typedef TaylorSeries<A> TS;
  typedef TaylorVariable<A> TV;


  class_<TS> taylor_series_class(python_name<R>("TaylorSeries").c_str());
  taylor_series_class.def("__init__", make_constructor(&make_taylor_series<R>));
  taylor_series_class.def( init< uint >());
  taylor_series_class.def( init< TV >());
  taylor_series_class.def("degree", &TS::degree);
  taylor_series_class.def("__getitem__", &get_item<A>);
  taylor_series_class.def("__setitem__",&set_item<A,double>);
  taylor_series_class.def("__setitem__",&set_item<A,R>);
  taylor_series_class.def("__setitem__",&set_item<A,A>);
  taylor_series_class.def("__neg__", &__neg__<TS,TS>);
  taylor_series_class.def("__add__", &__add__<TS,TS,TS>);
  taylor_series_class.def("__add__", &__add__<TS,TS,double>);
  //taylor_series_class.def("__add__", &__add__<TS,TS,R>);
  taylor_series_class.def("__add__", &__add__<TS,TS,A>);
  //taylor_series_class.def("__radd__", &__radd__<TS,TS,double>);
  //taylor_series_class.def("__radd__", &__radd__<TS,TS,R>);
  taylor_series_class.def("__radd__", &__radd__<TS,TS,A>);
  taylor_series_class.def("__sub__", &__sub__<TS,TS,TS>);
  taylor_series_class.def("__sub__", &__sub__<TS,TS,double>);
  ///taylor_series_class.def("__sub__", &__sub__<TS,TS,R>);
  taylor_series_class.def("__sub__", &__sub__<TS,TS,A>);
  //taylor_series_class.def("__rsub__", &__rsub__<TS,TS,double>);
  //taylor_series_class.def("__rsub__", &__rsub__<TS,TS,R>);
  taylor_series_class.def("__rsub__", &__rsub__<TS,TS,A>);
  taylor_series_class.def("__mul__", &__mul__<TS,TS,TS>);
  //taylor_series_class.def("__mul__", &__mul__<TS,TS,double>);
  //taylor_series_class.def("__mul__", &__mul__<TS,TS,R>);
  taylor_series_class.def("__mul__", &__mul__<TS,TS,A>);
  //taylor_series_class.def("__rmul__", &__rmul__<TS,TS,double>);
  //taylor_series_class.def("__rmul__", &__rmul__<TS,TS,R>);
  taylor_series_class.def("__rmul__", &__rmul__<TS,TS,A>);
  taylor_series_class.def("__div__", &__div__<TS,TS,TS>);
  //taylor_series_class.def("__div__", &__div__<TS,TS,double>);
  //taylor_series_class.def("__div__", &__div__<TS,TS,R>);
  taylor_series_class.def("__div__", &__div__<TS,TS,A>);
  //taylor_series_class.def("__rdiv__", &__rdiv__<TS,TS,double>);
  //taylor_series_class.def("__rdiv__", &__rdiv__<TS,TS,R>);
  taylor_series_class.def("__rdiv__", &__rdiv__<TS,TS,A>);
  taylor_series_class.def("__pow__", &__pow__<TS,TS,int>);
  taylor_series_class.def(self_ns::str(self));
  
  def("constant",(TS(*)(smoothness_type, const A&))&TS::constant);
  def("constant",(TS(*)(smoothness_type, const double&))&TS::constant);
  def("variable",(TS(*)(smoothness_type, const A&))&TS::variable);
  def("variable",(TS(*)(smoothness_type, const double&))&TS::variable);

  def("compose",(TS(*)(const TS&,const TS&))&compose);
  def("inverse",(TS(*)(const TS&,const A&))&inverse);


  def("max",(TS(*)(const TS&,const TS&))&max<A>);
  def("min",(TS(*)(const TS&,const TS&))&min<A>);
  def("abs",(TS(*)(const TS&))&abs<A>);
  def("pow",(TS(*)(const TS&, const uint&))&pow<A>);

  def("rec", (TS(*)(const TS&))&Ariadne::rec<A>);

  def("sqrt", (TS(*)(const TS&))&Ariadne::sqrt<A>);
  def("exp", (TS(*)(const TS&))&::exp<A>);
  def("log", (TS(*)(const TS&))&Ariadne::log<A>);
  def("sin", (TS(*)(const TS&))&Ariadne::sin<A>);
  def("cos", (TS(*)(const TS&))&Ariadne::cos<A>);
  def("tan", (TS(*)(const TS&))&Ariadne::tan<A>);
  def("asin", (TS(*)(const TS&))&Ariadne::asin<A>);
  def("acos", (TS(*)(const TS&))&Ariadne::acos<A>);
  def("atan", (TS(*)(const TS&))&Ariadne::atan<A>);

}

/*

template<>
void export_taylor_series<Rational>()
{
  typedef Rational Q;
  typedef TaylorSeries<Q> TS;


  class_<TS> taylor_series_class(python_name<Q>("TaylorSeries").c_str());
  taylor_series_class.def( init< uint >());
  taylor_series_class.def("__getitem__", &taylor_series_get_item<Q>);
  taylor_series_class.def("__setitem__",&taylor_series_set_item<Q,double>);
  taylor_series_class.def("__setitem__",&taylor_series_set_item<Q,Q>);
  taylor_series_class.def("__neg__", &__neg__<TS,TS>);
  taylor_series_class.def("__add__", &add<TS,TS,TS>);
  taylor_series_class.def("__add__", &add<TS,TS,double>);
  taylor_series_class.def("__add__", &add<TS,TS,Q>);
  taylor_series_class.def("__radd__", &radd<TS,TS,double>);
  taylor_series_class.def("__radd__", &radd<TS,TS,Q>);
  taylor_series_class.def("__sub__", &sub<TS,TS,TS>);
  taylor_series_class.def("__sub__", &sub<TS,TS,double>);
  taylor_series_class.def("__sub__", &sub<TS,TS,Q>);
  taylor_series_class.def("__rsub__", &rsub<TS,TS,double>);
  taylor_series_class.def("__rsub__", &rsub<TS,TS,Q>);
  taylor_series_class.def("__mul__", &mul<TS,TS,TS>);
  taylor_series_class.def("__mul__", &mul<TS,TS,double>);
  taylor_series_class.def("__mul__", &mul<TS,TS,Q>);
  taylor_series_class.def("__rmul__", &rmul<TS,TS,double>);
  taylor_series_class.def("__rmul__", &rmul<TS,TS,Q>);
  taylor_series_class.def("__div__", &div<TS,TS,TS>);
  taylor_series_class.def("__div__", &div<TS,TS,double>);
  taylor_series_class.def("__div__", &div<TS,TS,Q>);
  taylor_series_class.def("__rdiv__", &rdiv<TS,TS,double>);
  taylor_series_class.def("__rdiv__", &rdiv<TS,TS,Q>);
  taylor_series_class.def("__pow__", &pow<TS,TS,int>);
  taylor_series_class.def(self_ns::str(self));
  
  def("compose",(TS(*)(const TS&,const TS&))&compose);
  def("inverse",(TS(*)(const TS&,const Q&))&inverse);

  def("max",(TS(*)(const TS&,const TS&))&max);
  def("min",(TS(*)(const TS&,const TS&))&min);
  def("abs",(TS(*)(const TS&))&abs);
  def("pow",(TS(*)(const TS&, const uint&))&pow);

}

*/


 //template void export_taylor_series<Rational>();
template void export_taylor_series<FloatPy>();
