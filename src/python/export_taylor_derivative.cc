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
using namespace Ariadne::Function;
using namespace Ariadne::Python;


template<class X>
TaylorDerivative<X>*
make_taylor_derivative(const uint& rs, const uint& as, const uint& d, const boost::python::object& obj) 
{
  TaylorDerivative<X>* result=new TaylorDerivative<X>(rs,as,d);
  read_array(result->data(),obj);
  assert(result->data().size()==compute_polynomial_data_size(rs,as,d));
  return result;
}

template<class R1, class R2> inline 
void taylor_derivative_set_variable(TaylorDerivative<R1>& td, const size_type& i, TaylorVariable<R2> x) {
  assert(i==td.result_size()); 
  assert(x.argument_size()==td.argument_size()); 
  assert(x.degree()<=td.degree()); 
  td[i]=x;
}

template<class R> inline 
TaylorVariable<R> taylor_derivative_get_variable(const TaylorDerivative<R>& td, const size_type& i) {
  assert(i==td.result_size()); 
  return td[i];
}

template<class R1, class R2> inline 
void taylor_derivative_set_item(TaylorDerivative<R1>& td, const size_type& i, const MultiIndex& j, R2 x) {
  assert(i==td.result_size()); 
  assert(j.number_of_variables()==td.argument_size()); 
  assert(j.degree()<=td.degree()); 
  td.set(i,j,x);
}

template<class R> inline 
R taylor_derivative_get_item(const TaylorDerivative<R>& td, const size_type& i, const MultiIndex& j) {
  assert(i==td.result_size()); 
  assert(j.number_of_variables()==td.argument_size()); 
  assert(j.degree()<=td.degree()); 
  return td.get(i,j);
  ;
}

template<class R>
void export_taylor_derivative()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef TaylorVariable<A> TV;
  typedef TaylorDerivative<A> TD;

  class_<TD> taylor_derivative_class(python_name<R>("TaylorDerivative").c_str());
  taylor_derivative_class.def("__init__", make_constructor(&make_taylor_derivative<R>) );
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

}


template<>
void export_taylor_derivative<Rational>()
{
  typedef Rational Q;
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


}




template void export_taylor_derivative<Rational>();
template void export_taylor_derivative<FloatPy>();
