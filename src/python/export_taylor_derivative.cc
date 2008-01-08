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
 

#include "python/float.h"

#include "macros/assert.h"
#include "numeric/rational.h"
#include "numeric/interval.h"
#include "function/multi_index.h"
#include "function/taylor_derivative.h"

#include <boost/python.hpp>
#include "python/utilities.h"
#include "python/read_array.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Python;


template<class X>
TaylorDerivative<X>*
make_taylor_derivative(const uint& rs, const uint& as, const uint& d, const boost::python::object& obj) 
{
  array<X> data;
  read_array(data,obj);
  ARIADNE_ASSERT(data.size()==compute_polynomial_data_size(rs,as,d));
  TaylorDerivative<X>* result=new TaylorDerivative<X>(rs,as,d,data.begin());
  return result;
}

template<class X1, class X2> inline 
void taylor_derivative_set_variable(TaylorDerivative<X1>& td, const size_type& i, TaylorVariable<X2> x) {
  ARIADNE_ASSERT(i<td.result_size()); 
  ARIADNE_ASSERT(x.argument_size()==td.argument_size()); 
  ARIADNE_ASSERT(x.degree()<=td.degree()); 
  td[i]=x;
}

template<class X> inline 
TaylorVariable<X> taylor_derivative_get_variable(const TaylorDerivative<X>& td, const size_type& i) {
  ARIADNE_ASSERT(i<td.result_size()); 
  return td[i];
}

template<class X, class XX> inline 
void taylor_derivative_set_item(TaylorDerivative<X>& td, const size_type& i, const MultiIndex& j, const XX& x) {
  ARIADNE_ASSERT(i<td.result_size()); 
  ARIADNE_ASSERT(j.number_of_variables()==td.argument_size()); 
  ARIADNE_ASSERT(j.degree()<=td.degree()); 
  td.set(i,j,x);
}

template<class X> inline 
X taylor_derivative_get_item(const TaylorDerivative<X>& td, const size_type& i, const MultiIndex& j) {
  ARIADNE_ASSERT(i==td.result_size()); 
  ARIADNE_ASSERT(j.number_of_variables()==td.argument_size()); 
  ARIADNE_ASSERT(j.degree()<=td.degree()); 
  return td.get(i,j);
  ;
}

template<class R>
void export_taylor_derivative()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef Vector<A> Vec;
  typedef TaylorVariable<A> TV;
  typedef TaylorDerivative<A> TD;

  class_<TD> taylor_derivative_class(python_name<R>("TaylorDerivative").c_str());
  taylor_derivative_class.def("__init__", make_constructor(&make_taylor_derivative<A>) );
  taylor_derivative_class.def( init< uint, uint, uint >());
  taylor_derivative_class.def("__getitem__", &taylor_derivative_get_item<A>);
  taylor_derivative_class.def("__setitem__",&taylor_derivative_set_item<A,double>);
  taylor_derivative_class.def("__setitem__",&taylor_derivative_set_item<A,R>);
  taylor_derivative_class.def("__setitem__",&taylor_derivative_set_item<A,A>);
  taylor_derivative_class.def("__getitem__", &taylor_derivative_get_variable<A>);
  taylor_derivative_class.def("__setitem__",&taylor_derivative_set_variable<A,A>);
  taylor_derivative_class.def("__neg__", &Python::neg<TD,TD>);
  taylor_derivative_class.def("__add__", &Python::add<TD,TD,TD>);
  taylor_derivative_class.def("__sub__", &Python::sub<TD,TD,TD>);
  taylor_derivative_class.def(self_ns::str(self));
  
  def("compose",(TV(*)(const TV&,const TD&))&Function::compose);
  def("compose",(TD(*)(const TD&,const TD&))&Function::compose);
  def("inverse",(TD(*)(const TD&,const Vec&))&Function::inverse);
  def("implicit",(TD(*)(const TD&,const Vec&))&Function::implicit);

}


template<>
void export_taylor_derivative<Rational>()
{
  typedef Rational Q;
  typedef Vector<Q> Vec;
  typedef TaylorVariable<Q> TV;
  typedef TaylorDerivative<Q> TD;


  class_<TD> taylor_derivative_class(python_name<Q>("TaylorDerivative").c_str());
  taylor_derivative_class.def( init< uint,uint, uint >());
  taylor_derivative_class.def("__getitem__", &taylor_derivative_get_item<Q>);
  taylor_derivative_class.def("__setitem__",&taylor_derivative_set_item<Q,double>);
  taylor_derivative_class.def("__setitem__",&taylor_derivative_set_item<Q,Q>);
  taylor_derivative_class.def("__getitem__", &taylor_derivative_get_variable<Q>);
  taylor_derivative_class.def("__setitem__",&taylor_derivative_set_variable<Q,Q>);
  taylor_derivative_class.def("__neg__", &Python::neg<TD,TD>);
  taylor_derivative_class.def("__add__", &Python::add<TD,TD,TD>);
  taylor_derivative_class.def("__sub__", &Python::sub<TD,TD,TD>);
  taylor_derivative_class.def(self_ns::str(self));
  
  def("compose",(TV(*)(const TV&,const TD&))&Function::compose);
  def("compose",(TD(*)(const TD&,const TD&))&Function::compose);
  def("inverse",(TD(*)(const TD&,const Vec&))&Function::inverse);
  def("implicit",(TD(*)(const TD&,const Vec&))&Function::implicit);

  
}




template void export_taylor_derivative<Rational>();
template void export_taylor_derivative<FloatPy>();
