/***************************************************************************
 *            python/export_affine_variable.cc
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
#include "function/affine_variable.h"

#include <boost/python.hpp>
#include "python/python_utilities.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace Ariadne::Python;


template<class R1, class R2> inline 
void affine_variable_set_item(AffineVariable<R1>& td, const MultiIndex& i, R2 x) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  td[i]=x;
}

template<class R> inline 
R affine_variable_get_item(const AffineVariable<R>& td, const MultiIndex& i) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  return td[i];
}

template<class R>
void export_affine_variable()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef ScalarDerivative<A> SD;
  typedef AffineVariable<A> TD;


  class_<TD> affine_variable_class(python_name<R>("AffineVariable").c_str());
  affine_variable_class.def( init< uint, uint >());
  affine_variable_class.def( init< uint, uint, double >());
  affine_variable_class.def( init< uint, uint, uint, double >());
  affine_variable_class.def( init< uint, uint, A >());
  affine_variable_class.def( init< uint, uint, uint, A >());
  affine_variable_class.def("__getitem__", &affine_variable_get_item<A>);
  affine_variable_class.def("__setitem__",&affine_variable_set_item<A,double>);
  affine_variable_class.def("__setitem__",&affine_variable_set_item<A,R>);
  affine_variable_class.def("__setitem__",&affine_variable_set_item<A,A>);
  affine_variable_class.def("__neg__", &Python::neg<TD,TD>);
  affine_variable_class.def("__add__", &add<TD,TD,TD>);
  affine_variable_class.def("__add__", &add<TD,TD,double>);
  affine_variable_class.def("__add__", &add<TD,TD,R>);
  affine_variable_class.def("__add__", &add<TD,TD,A>);
  affine_variable_class.def("__radd__", &radd<TD,TD,double>);
  affine_variable_class.def("__radd__", &radd<TD,TD,R>);
  affine_variable_class.def("__radd__", &radd<TD,TD,A>);
  affine_variable_class.def("__sub__", &sub<TD,TD,TD>);
  affine_variable_class.def("__sub__", &sub<TD,TD,double>);
  affine_variable_class.def("__sub__", &sub<TD,TD,R>);
  affine_variable_class.def("__sub__", &sub<TD,TD,A>);
  affine_variable_class.def("__rsub__", &rsub<TD,TD,double>);
  affine_variable_class.def("__rsub__", &rsub<TD,TD,R>);
  affine_variable_class.def("__rsub__", &rsub<TD,TD,A>);
  affine_variable_class.def("__mul__", &mul<TD,TD,TD>);
  affine_variable_class.def("__mul__", &mul<TD,TD,double>);
  affine_variable_class.def("__mul__", &mul<TD,TD,R>);
  affine_variable_class.def("__mul__", &mul<TD,TD,A>);
  affine_variable_class.def("__rmul__", &rmul<TD,TD,double>);
  affine_variable_class.def("__rmul__", &rmul<TD,TD,R>);
  affine_variable_class.def("__rmul__", &rmul<TD,TD,A>);
  affine_variable_class.def("__div__", &div<TD,TD,TD>);
  affine_variable_class.def("__div__", &div<TD,TD,double>);
  affine_variable_class.def("__div__", &div<TD,TD,R>);
  affine_variable_class.def("__div__", &div<TD,TD,A>);
  affine_variable_class.def("__rdiv__", &rdiv<TD,TD,double>);
  affine_variable_class.def("__rdiv__", &rdiv<TD,TD,R>);
  affine_variable_class.def("__rdiv__", &rdiv<TD,TD,A>);
  affine_variable_class.def("__pow__", &pow<TD,TD,int>);
  affine_variable_class.def(self_ns::str(self));

  
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
void export_affine_variable<Rational>()
{
  typedef Rational Q;
  typedef ScalarDerivative<Q> SD;
  typedef AffineVariable<Q> TD;


  class_<TD>  affine_variable_class(python_name<Q>("AffineVariable").c_str());
  affine_variable_class.def( init< uint, uint >());
  affine_variable_class.def( init< uint, uint, double >());
  affine_variable_class.def( init< uint, uint, Q >());
  affine_variable_class.def( init< uint, uint, uint, double >());
  affine_variable_class.def( init< uint, uint, uint, Q >());
  affine_variable_class.def("__getitem__", &affine_variable_get_item<Q>);
  affine_variable_class.def("__setitem__",&affine_variable_set_item<Q,double>);
  affine_variable_class.def("__setitem__",&affine_variable_set_item<Q,Q>);
  affine_variable_class.def("__neg__", &Python::neg<TD,TD>);
  affine_variable_class.def("__add__", &add<TD,TD,TD>);
  affine_variable_class.def("__add__", &add<TD,TD,double>);
  affine_variable_class.def("__add__", &add<TD,TD,Q>);
  affine_variable_class.def("__radd__", &radd<TD,TD,double>);
  affine_variable_class.def("__radd__", &radd<TD,TD,Q>);
  affine_variable_class.def("__sub__", &sub<TD,TD,TD>);
  affine_variable_class.def("__sub__", &sub<TD,TD,double>);
  affine_variable_class.def("__sub__", &sub<TD,TD,Q>);
  affine_variable_class.def("__rsub__", &rsub<TD,TD,double>);
  affine_variable_class.def("__rsub__", &rsub<TD,TD,Q>);
  affine_variable_class.def("__mul__", &mul<TD,TD,TD>);
  affine_variable_class.def("__mul__", &mul<TD,TD,double>);
  affine_variable_class.def("__mul__", &mul<TD,TD,Q>);
  affine_variable_class.def("__rmul__", &rmul<TD,TD,double>);
  affine_variable_class.def("__rmul__", &rmul<TD,TD,Q>);
  affine_variable_class.def("__div__", &div<TD,TD,TD>);
  affine_variable_class.def("__div__", &div<TD,TD,double>);
  affine_variable_class.def("__div__", &div<TD,TD,Q>);
  affine_variable_class.def("__rdiv__", &rdiv<TD,TD,double>);
  affine_variable_class.def("__rdiv__", &rdiv<TD,TD,Q>);
  affine_variable_class.def("__pow__", &pow<TD,TD,int>);
  affine_variable_class.def(self_ns::str(self));

  
  def("compose",(TD(*)(const SD&,const TD&))&Function::compose);

  def("max",(TD(*)(const TD&,const TD&))&Function::max);
  def("min",(TD(*)(const TD&,const TD&))&Function::min);
  def("abs",(TD(*)(const TD&))&Function::abs);
  def("pow",(TD(*)(const TD&, int))&Function::pow);

}




template void export_affine_variable<Rational>();
template void export_affine_variable<Float>();
