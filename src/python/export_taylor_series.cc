/***************************************************************************
 *            python/export_taylor_series`.cc
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
#include "function/taylor_series.h"
#include "function/taylor_variable.h"

#include <boost/python.hpp>
#include "python/utilities.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace Ariadne::Python;


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

template<class R>
void export_taylor_series()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef TaylorSeries<A> TS;
  typedef TaylorVariable<A> TV;


  class_<TS> taylor_series_class(python_name<R>("TaylorSeries").c_str());
  taylor_series_class.def( init< uint >());
  taylor_series_class.def( init< TV >());
  taylor_series_class.def("degree", &TS::degree);
  taylor_series_class.def("__getitem__", &get_item<A>);
  taylor_series_class.def("__setitem__",&set_item<A,double>);
  taylor_series_class.def("__setitem__",&set_item<A,R>);
  taylor_series_class.def("__setitem__",&set_item<A,A>);
  taylor_series_class.def("__neg__", &Python::neg<TS,TS>);
  taylor_series_class.def("__add__", &Python::add<TS,TS,TS>);
  //taylor_series_class.def("__add__", &Python::add<TS,TS,double>);
  //taylor_series_class.def("__add__", &Python::add<TS,TS,R>);
  taylor_series_class.def("__add__", &Python::add<TS,TS,A>);
  //taylor_series_class.def("__radd__", &Python::radd<TS,TS,double>);
  //taylor_series_class.def("__radd__", &Python::radd<TS,TS,R>);
  taylor_series_class.def("__radd__", &Python::radd<TS,TS,A>);
  taylor_series_class.def("__sub__", &Python::sub<TS,TS,TS>);
  //taylor_series_class.def("__sub__", &Python::sub<TS,TS,double>);
  ///taylor_series_class.def("__sub__", &Python::sub<TS,TS,R>);
  taylor_series_class.def("__sub__", &Python::sub<TS,TS,A>);
  //taylor_series_class.def("__rsub__", &Python::rsub<TS,TS,double>);
  //taylor_series_class.def("__rsub__", &Python::rsub<TS,TS,R>);
  taylor_series_class.def("__rsub__", &Python::rsub<TS,TS,A>);
  taylor_series_class.def("__mul__", &Python::mul<TS,TS,TS>);
  //taylor_series_class.def("__mul__", &Python::mul<TS,TS,double>);
  //taylor_series_class.def("__mul__", &Python::mul<TS,TS,R>);
  taylor_series_class.def("__mul__", &Python::mul<TS,TS,A>);
  //taylor_series_class.def("__rmul__", &Python::rmul<TS,TS,double>);
  //taylor_series_class.def("__rmul__", &Python::rmul<TS,TS,R>);
  taylor_series_class.def("__rmul__", &Python::rmul<TS,TS,A>);
  taylor_series_class.def("__div__", &Python::div<TS,TS,TS>);
  //taylor_series_class.def("__div__", &Python::div<TS,TS,double>);
  //taylor_series_class.def("__div__", &Python::div<TS,TS,R>);
  taylor_series_class.def("__div__", &Python::div<TS,TS,A>);
  //taylor_series_class.def("__rdiv__", &Python::rdiv<TS,TS,double>);
  //taylor_series_class.def("__rdiv__", &Python::rdiv<TS,TS,R>);
  taylor_series_class.def("__rdiv__", &Python::rdiv<TS,TS,A>);
  taylor_series_class.def("__pow__", &Python::pow<TS,TS,int>);
  taylor_series_class.def(self_ns::str(self));
  
  def("constant",(TS(*)(smoothness_type, const A&))&TS::constant);
  def("constant",(TS(*)(smoothness_type, const double&))&TS::constant);
  def("variable",(TS(*)(smoothness_type, const A&))&TS::variable);
  def("variable",(TS(*)(smoothness_type, const double&))&TS::variable);

  def("compose",(TS(*)(const TS&,const TS&))&Function::compose);
  def("inverse",(TS(*)(const TS&,const A&))&Function::inverse);

  def("max",(TS(*)(const TS&,const TS&))&Function::max);
  def("min",(TS(*)(const TS&,const TS&))&Function::min);
  def("abs",(TS(*)(const TS&))&Function::abs);
  def("pow",(TS(*)(const TS&, const uint&))&Function::pow);

  def("rec", (TS(*)(const TS&))&Function::rec);

  def("sqrt", (TS(*)(const TS&))&Function::sqrt);
  def("exp", (TS(*)(const TS&))&Function::exp);
  def("log", (TS(*)(const TS&))&Function::log);
  def("sin", (TS(*)(const TS&))&Function::sin);
  def("cos", (TS(*)(const TS&))&Function::cos);
  def("tan", (TS(*)(const TS&))&Function::tan);
  def("asin", (TS(*)(const TS&))&Function::asin);
  def("acos", (TS(*)(const TS&))&Function::acos);
  def("atan", (TS(*)(const TS&))&Function::atan);

  /*
  def("rec_series", (TS(*)(smoothness_type, const A&))&TS::rec);
  def("pow_series",(TS(*)(smoothness_type, const A&, const uint&))&TS::pow);
  def("sqrt_series", (TS(*)(smoothness_type, const A&))&TS::sqrt);
  def("exp_series", (TS(*)(smoothness_type, const A&))&TS::exp);
  def("log_series", (TS(*)(smoothness_type, const A&))&TS::log);
  def("sin_series", (TS(*)(smoothness_type, const A&))&TS::sin);
  def("cos_series", (TS(*)(smoothness_type, const A&))&TS::cos);
  def("tan_series", (TS(*)(smoothness_type, const A&))&TS::tan);
  def("asin_series", (TS(*)(smoothness_type, const A&))&TS::asin);
  def("acos_series", (TS(*)(smoothness_type, const A&))&TS::acos);
  def("atan_series", (TS(*)(smoothness_type, const A&))&TS::atan);
  */

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
  taylor_series_class.def("__neg__", &Python::neg<TS,TS>);
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
  
  def("compose",(TS(*)(const TS&,const TS&))&Function::compose);
  def("inverse",(TS(*)(const TS&,const Q&))&Function::inverse);

  def("max",(TS(*)(const TS&,const TS&))&Function::max);
  def("min",(TS(*)(const TS&,const TS&))&Function::min);
  def("abs",(TS(*)(const TS&))&Function::abs);
  def("pow",(TS(*)(const TS&, const uint&))&Function::pow);

}

*/


 //template void export_taylor_series<Rational>();
template void export_taylor_series<FloatPy>();
