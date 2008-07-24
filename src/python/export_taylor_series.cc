/***************************************************************************
 *            python/export_power_series.cc
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
#include "differentiation/power_series.h"
#include "differentiation/differential.h"
#include "python/read_array.h"
using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include "python/utilities.h"
using namespace boost::python;


template<class X>
PowerSeries<X>*
make_power_series(const boost::python::object& obj) 
{
  PowerSeries<X>* ts=new PowerSeries<X>;
  read_array(ts->data(),obj);
  return ts;
}

template<class R1, class R2> inline 
void set_item(PowerSeries<R1>& sd, uint i, R2 x) {
  assert(i<=sd.degree()); 
  sd[i]=x;
}

template<class R> inline 
R get_item(const PowerSeries<R>& sd, uint i) {
  assert(i<=sd.degree()); 
  return sd[i];
}

template<class X> inline
PowerSeries<X> exp(const PowerSeries<X>& ts) {
  return Ariadne::exp(ts);
}


template<class R>
void export_power_series()
{
  typedef typename traits<R>::arithmetic_type A;
  typedef PowerSeries<A> TS;
  typedef Differential<A> TV;


  class_<TS> power_series_class(python_name<R>("PowerSeries").c_str());
  power_series_class.def("__init__", make_constructor(&make_power_series<R>));
  power_series_class.def( init< uint >());
  power_series_class.def("degree", &TS::degree);
  power_series_class.def("__getitem__", &get_item<A>);
  power_series_class.def("__setitem__",&set_item<A,double>);
  power_series_class.def("__setitem__",&set_item<A,R>);
  power_series_class.def("__setitem__",&set_item<A,A>);
  power_series_class.def("__neg__", &__neg__<TS,TS>);
  power_series_class.def("__add__", &__add__<TS,TS,TS>);
  power_series_class.def("__add__", &__add__<TS,TS,double>);
  //power_series_class.def("__add__", &__add__<TS,TS,R>);
  power_series_class.def("__add__", &__add__<TS,TS,A>);
  //power_series_class.def("__radd__", &__radd__<TS,TS,double>);
  //power_series_class.def("__radd__", &__radd__<TS,TS,R>);
  power_series_class.def("__radd__", &__radd__<TS,TS,A>);
  power_series_class.def("__sub__", &__sub__<TS,TS,TS>);
  power_series_class.def("__sub__", &__sub__<TS,TS,double>);
  ///power_series_class.def("__sub__", &__sub__<TS,TS,R>);
  power_series_class.def("__sub__", &__sub__<TS,TS,A>);
  //power_series_class.def("__rsub__", &__rsub__<TS,TS,double>);
  //power_series_class.def("__rsub__", &__rsub__<TS,TS,R>);
  power_series_class.def("__rsub__", &__rsub__<TS,TS,A>);
  power_series_class.def("__mul__", &__mul__<TS,TS,TS>);
  //power_series_class.def("__mul__", &__mul__<TS,TS,double>);
  //power_series_class.def("__mul__", &__mul__<TS,TS,R>);
  power_series_class.def("__mul__", &__mul__<TS,TS,A>);
  //power_series_class.def("__rmul__", &__rmul__<TS,TS,double>);
  //power_series_class.def("__rmul__", &__rmul__<TS,TS,R>);
  power_series_class.def("__rmul__", &__rmul__<TS,TS,A>);
  power_series_class.def("__div__", &__div__<TS,TS,TS>);
  //power_series_class.def("__div__", &__div__<TS,TS,double>);
  //power_series_class.def("__div__", &__div__<TS,TS,R>);
  power_series_class.def("__div__", &__div__<TS,TS,A>);
  //power_series_class.def("__rdiv__", &__rdiv__<TS,TS,double>);
  //power_series_class.def("__rdiv__", &__rdiv__<TS,TS,R>);
  power_series_class.def("__rdiv__", &__rdiv__<TS,TS,A>);
  power_series_class.def("__pow__", &__pow__<TS,TS,int>);
  power_series_class.def(self_ns::str(self));
  
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
void export_power_series<Rational>()
{
  typedef Rational Q;
  typedef PowerSeries<Q> TS;


  class_<TS> power_series_class(python_name<Q>("PowerSeries").c_str());
  power_series_class.def( init< uint >());
  power_series_class.def("__getitem__", &power_series_get_item<Q>);
  power_series_class.def("__setitem__",&power_series_set_item<Q,double>);
  power_series_class.def("__setitem__",&power_series_set_item<Q,Q>);
  power_series_class.def("__neg__", &__neg__<TS,TS>);
  power_series_class.def("__add__", &add<TS,TS,TS>);
  power_series_class.def("__add__", &add<TS,TS,double>);
  power_series_class.def("__add__", &add<TS,TS,Q>);
  power_series_class.def("__radd__", &radd<TS,TS,double>);
  power_series_class.def("__radd__", &radd<TS,TS,Q>);
  power_series_class.def("__sub__", &sub<TS,TS,TS>);
  power_series_class.def("__sub__", &sub<TS,TS,double>);
  power_series_class.def("__sub__", &sub<TS,TS,Q>);
  power_series_class.def("__rsub__", &rsub<TS,TS,double>);
  power_series_class.def("__rsub__", &rsub<TS,TS,Q>);
  power_series_class.def("__mul__", &mul<TS,TS,TS>);
  power_series_class.def("__mul__", &mul<TS,TS,double>);
  power_series_class.def("__mul__", &mul<TS,TS,Q>);
  power_series_class.def("__rmul__", &rmul<TS,TS,double>);
  power_series_class.def("__rmul__", &rmul<TS,TS,Q>);
  power_series_class.def("__div__", &div<TS,TS,TS>);
  power_series_class.def("__div__", &div<TS,TS,double>);
  power_series_class.def("__div__", &div<TS,TS,Q>);
  power_series_class.def("__rdiv__", &rdiv<TS,TS,double>);
  power_series_class.def("__rdiv__", &rdiv<TS,TS,Q>);
  power_series_class.def("__pow__", &pow<TS,TS,int>);
  power_series_class.def(self_ns::str(self));
  
  def("compose",(TS(*)(const TS&,const TS&))&compose);
  def("inverse",(TS(*)(const TS&,const Q&))&inverse);

  def("max",(TS(*)(const TS&,const TS&))&max);
  def("min",(TS(*)(const TS&,const TS&))&min);
  def("abs",(TS(*)(const TS&))&abs);
  def("pow",(TS(*)(const TS&, const uint&))&pow);

}

*/


 //template void export_power_series<Rational>();
template void export_power_series<FloatPy>();
