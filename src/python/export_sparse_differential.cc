/***************************************************************************
 *            python/export_sparse_differential.cc
 *
 *  Copyright  2008 Pieter Collins
 *
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
#include "differentiation/sparse_differential.h"
#include "differentiation/taylor_variable.h"

#include "python/utilities.h"
#include "python/read_array.h"
using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class X>
SparseDifferential<X>*
make_sparse_differential(const uint& as, const uint& d, const boost::python::object& obj) 
{
  array<X> ary; 
  read_array(ary,obj);
  SparseDifferential<X>* result=new SparseDifferential<X>(as,d,ary.begin());
  return result;
}

template<class X, class XX> inline 
void sparse_differential_set_item(SparseDifferential<X>& sd, const MultiIndex& i, const XX& x) {
  assert(i.number_of_variables()==sd.argument_size()); 
  assert(i.degree()<=sd.degree()); 
  sd[i]=X(x);
}

template<class X> inline 
X sparse_differential_get_item(const SparseDifferential<X>& sd, const MultiIndex& i) {
  assert(i.number_of_variables()==sd.argument_size()); 
  assert(i.degree()<=sd.degree()); 
  return sd[i];
}

template<class X>
SparseDifferential<X> exp(const SparseDifferential<X>& sd) {
  return Ariadne::exp(sd);
}


template<class X>
SparseDifferentialVector<X>*
make_sparse_differential_vector(const uint& rs, const uint& as, const uint& d, const boost::python::object& obj) 
{
  array<X> data;
  read_array(data,obj);
  ARIADNE_ASSERT(data.size()==rs*compute_polynomial_data_size(as,d));
  SparseDifferentialVector<X>* result=new SparseDifferentialVector<X>(rs,as,d,data.begin());
  return result;
}

template<class X1, class X2> inline 
void sparse_differential_vector_set_variable(SparseDifferentialVector<X1>& sdv, const size_type& i, SparseDifferential<X2> x) {
  ARIADNE_ASSERT(i<sdv.result_size()); 
  ARIADNE_ASSERT(x.argument_size()==sdv.argument_size()); 
  ARIADNE_ASSERT(x.degree()<=sdv.degree()); 
  sdv[i]=x;
}

template<class X> inline 
SparseDifferential<X> sparse_differential_vector_get_variable(const SparseDifferentialVector<X>& sdv, const size_type& i) {
  ARIADNE_ASSERT(i<sdv.result_size()); 
  return sdv[i];
}

template<class X, class XX> inline 
void sparse_differential_vector_set_item(SparseDifferentialVector<X>& sdv, const size_type& i, const MultiIndex& j, const XX& x) {
  ARIADNE_ASSERT(i<sdv.result_size()); 
  ARIADNE_ASSERT(j.number_of_variables()==sdv.argument_size()); 
  ARIADNE_ASSERT(j.degree()<=sdv.degree()); 
  sdv[i][j]=X(x);
}

template<class X> inline 
X sparse_differential_vector_get_item(const SparseDifferentialVector<X>& sdv, const size_type& i, const MultiIndex& j) {
  ARIADNE_ASSERT(i==sdv.result_size()); 
  ARIADNE_ASSERT(j.number_of_variables()==sdv.argument_size()); 
  ARIADNE_ASSERT(j.degree()<=sdv.degree()); 
  return sdv[i][j];
}

template<class X, class XX> inline 
void sparse_differential_vector_set_gradient(SparseDifferentialVector<X>& sdv, const size_type& i, const size_type& j, const XX& x) {
  ARIADNE_ASSERT(i<sdv.result_size()); 
  ARIADNE_ASSERT(j<=sdv.argument_size());
  ARIADNE_ASSERT(1u<=sdv.degree()); 
  sdv[i][j]=X(x);
}

template<class X> inline 
X sparse_differential_vector_get_gradient(const SparseDifferentialVector<X>& sdv, const size_type& i, const size_type& j) {
  ARIADNE_ASSERT(i<sdv.result_size()); 
  ARIADNE_ASSERT(j<=sdv.argument_size());
  ARIADNE_ASSERT(1u<=sdv.degree()); 
  return sdv[i][j];
}


template<class R>
void export_sparse_differential()
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef TaylorSeries<A> TS;
  typedef SparseDifferential<A> SD;


  class_<SD> sparse_differential_class(python_name<R>("SparseDifferential").c_str());
  sparse_differential_class.def("__init__", make_constructor(&make_sparse_differential<A>) );
  sparse_differential_class.def( init< uint, uint >());
  sparse_differential_class.def("__getitem__", &sparse_differential_get_item<A>);
  sparse_differential_class.def("__setitem__",&sparse_differential_set_item<A,double>);
  sparse_differential_class.def("__setitem__",&sparse_differential_set_item<A,R>);
  sparse_differential_class.def("__setitem__",&sparse_differential_set_item<A,A>);
  sparse_differential_class.def("__neg__", &__neg__<SD,SD>);
  sparse_differential_class.def("__add__", &__add__<SD,SD,SD>);
  sparse_differential_class.def("__add__", &__add__<SD,SD,double>);
  sparse_differential_class.def("__add__", &__add__<SD,SD,R>);
  sparse_differential_class.def("__add__", &__add__<SD,SD,A>);
  sparse_differential_class.def("__radd__", &__radd__<SD,SD,double>);
  sparse_differential_class.def("__radd__", &__radd__<SD,SD,R>);
  sparse_differential_class.def("__radd__", &__radd__<SD,SD,A>);
  sparse_differential_class.def("__sub__", &__sub__<SD,SD,SD>);
  sparse_differential_class.def("__sub__", &__sub__<SD,SD,double>);
  sparse_differential_class.def("__sub__", &__sub__<SD,SD,R>);
  sparse_differential_class.def("__sub__", &__sub__<SD,SD,A>);
  sparse_differential_class.def("__rsub__", &__rsub__<SD,SD,double>);
  sparse_differential_class.def("__rsub__", &__rsub__<SD,SD,R>);
  sparse_differential_class.def("__rsub__", &__rsub__<SD,SD,A>);
  sparse_differential_class.def("__mul__", &__mul__<SD,SD,SD>);
  sparse_differential_class.def("__mul__", &__mul__<SD,SD,double>);
  sparse_differential_class.def("__mul__", &__mul__<SD,SD,R>);
  sparse_differential_class.def("__mul__", &__mul__<SD,SD,A>);
  sparse_differential_class.def("__rmul__", &__rmul__<SD,SD,double>);
  sparse_differential_class.def("__rmul__", &__rmul__<SD,SD,R>);
  sparse_differential_class.def("__rmul__", &__rmul__<SD,SD,A>);
  sparse_differential_class.def("__div__", &__div__<SD,SD,SD>);
  sparse_differential_class.def("__div__", &__div__<SD,SD,double>);
  sparse_differential_class.def("__div__", &__div__<SD,SD,R>);
  sparse_differential_class.def("__div__", &__div__<SD,SD,A>);
  sparse_differential_class.def("__rdiv__", &__rdiv__<SD,SD,double>);
  //sparse_differential_class.def("__rdiv__", &__rdiv__<SD,SD,R>);
  sparse_differential_class.def("__rdiv__", &__rdiv__<SD,SD,A>);
  sparse_differential_class.def("__pow__", &__pow__<SD,SD,int>);
  sparse_differential_class.def(self_ns::str(self));
  
  def("sparse_constant",(SD(*)(size_type, smoothness_type, const double&))&SD::constant);
  def("sparse_variable",(SD(*)(size_type, smoothness_type, const double&, size_type))&SD::variable);

  def("sparse_constant",(SD(*)(size_type, smoothness_type, const A&))&SD::constant);
  def("sparse_variable",(SD(*)(size_type, smoothness_type, const A&, size_type))&SD::variable);

  def("compose",(SD(*)(const TS&,const SD&))&compose);

  //def("max",(SD(*)(const SD&,const SD&))&max<A>);
  //def("min",(SD(*)(const SD&,const SD&))&min<A>);
  //def("abs",(SD(*)(const SD&))&abs<A>);
  //def("neg",(SD(*)(const SD&))&neg<A>);
  def("rec",(SD(*)(const SD&))&rec<A>);
  def("pow",(SD(*)(const SD&, int))&pow<A>);

  def("sqrt", (SD(*)(const SD&))&sqrt<A>);
  def("exp", (SD(*)(const SD&))&::exp<A>);
  def("log", (SD(*)(const SD&))&log<A>);
  def("sin", (SD(*)(const SD&))&sin<A>);
  def("cos", (SD(*)(const SD&))&cos<A>);
  def("tan", (SD(*)(const SD&))&tan<A>);
  //def("asin", (SD(*)(const SD&))&asin<A>);
  //def("acos", (SD(*)(const SD&))&acos<A>);
  //def("atan", (SD(*)(const SD&))&atan<A>);

}


template<class R>
void export_sparse_differential_vector()
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef typename traits<R>::interval_type I;
  typedef Vector<R> RVec;
  typedef Vector<A> AVec;
  typedef Vector<I> IVec;
  typedef SparseDifferential<A> SD;
  typedef SparseDifferentialVector<A> SDV;

  class_<SDV> sparse_differential_vector_class(python_name<R>("SparseDifferentialVector").c_str());
  sparse_differential_vector_class.def("__init__", make_constructor(&make_sparse_differential_vector<A>) );
  sparse_differential_vector_class.def( init< uint, uint, uint >());
  sparse_differential_vector_class.def("__getitem__", &sparse_differential_vector_get_item<A>);
  sparse_differential_vector_class.def("__setitem__",&sparse_differential_vector_set_item<A,double>);
  sparse_differential_vector_class.def("__setitem__",&sparse_differential_vector_set_item<A,R>);
  sparse_differential_vector_class.def("__setitem__",&sparse_differential_vector_set_item<A,A>);
  sparse_differential_vector_class.def("__getitem__", &sparse_differential_vector_get_variable<A>);
  sparse_differential_vector_class.def("__setitem__",&sparse_differential_vector_set_variable<A,A>);
  sparse_differential_vector_class.def("__setitem__",&sparse_differential_vector_set_gradient<A,double>);
  //sparse_differential_vector_class.def("__neg__", &__neg__<SDV,SDV>);
  sparse_differential_vector_class.def("__add__", &__add__<SDV,SDV,SDV>);
  sparse_differential_vector_class.def("__sub__", &__sub__<SDV,SDV,SDV>);
  sparse_differential_vector_class.def("__sub__", &__sub__<SDV,SDV,AVec>);
  //sparse_differential_vector_class.def("__mul__", &__mul__<SDV,SDV,A>);
  sparse_differential_vector_class.def("value", &SDV::value);
  sparse_differential_vector_class.def("jacobian", &SDV::jacobian);
  sparse_differential_vector_class.def(self_ns::str(self));
  
  def("sparse_variable",(SDV(*)(const RVec&,smoothness_type))&SDV::variable);
  def("sparse_variable",(SDV(*)(const IVec&,smoothness_type))&SDV::variable);

  def("evaluate",(IVec(*)(const SDV&,const IVec&))&evaluate);
  //def("evaluate",(SD(*)(const SD&,const SDV&))&evaluate);
  //def("evaluate",(SDV(*)(const SDV&,const SDV&))&evaluate);
  def("compose",(SD(*)(const SD&,const SDV&))&compose);
  def("compose",(SDV(*)(const SDV&,const SDV&))&compose);
  def("translate",(SDV(*)(const SDV&,const AVec&))&translate);
  def("inverse",(SDV(*)(const SDV&))&inverse);
  def("implicit",(SDV(*)(const SDV&))&implicit);
  def("flow",(SDV(*)(const SDV&,const AVec&,ushort,ushort))&flow);

}



template void export_sparse_differential<FloatPy>();
template void export_sparse_differential_vector<FloatPy>();
