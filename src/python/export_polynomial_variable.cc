/***************************************************************************
 *            python/export_polynomial_variable.cc
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
#include "function/polynomial_variable.h"

#include <boost/python.hpp>
#include "python/utilities.h"
#include "python/read_array.h"

using namespace boost::python;
using namespace Ariadne;

using namespace Ariadne::Function;
using namespace Ariadne::Python;

template<class X>
PolynomialVariable<X>*
make_polynomial_variable(const uint& as, const uint& d, const boost::python::object& obj) 
{
  PolynomialVariable<X>* result=new PolynomialVariable<X>(as,d);
  read_array(result->data(),obj);
  assert(result->data().size()==compute_polynomial_data_size(as,d));
  return result;
}


template<class R1, class R2> inline 
void polynomial_variable_set(PolynomialVariable<R1>& td, const uint& i, R2 x) {
  td.data()[i]=x;
}

template<class R1, class R2> inline 
void polynomial_variable_set_item(PolynomialVariable<R1>& td, const MultiIndex& i, R2 x) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  td[i]=x;
}

template<class R> inline 
R polynomial_variable_get_item(const PolynomialVariable<R>& td, const MultiIndex& i) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  return td[i];
}

template<class R>
void export_polynomial_variable()
{
  typedef typename traits<R>::arithmetic_type A;
  typedef PolynomialVariable<A> PV;


  class_<PV> polynomial_variable_class(python_name<R>("PolynomialVariable").c_str());
  polynomial_variable_class.def("__init__", make_constructor(&make_polynomial_variable<A>) );
  polynomial_variable_class.def( init< uint, uint >());
  polynomial_variable_class.def("__getitem__", &polynomial_variable_get_item<A>);
  polynomial_variable_class.def("__setitem__",&polynomial_variable_set<A,double>);
  polynomial_variable_class.def("__setitem__",&polynomial_variable_set_item<A,R>);
  polynomial_variable_class.def("__setitem__",&polynomial_variable_set_item<A,A>);
  polynomial_variable_class.def("__neg__", &Python::neg<PV,PV>);
  polynomial_variable_class.def("__add__", &Python::add<PV,PV,PV>);
  polynomial_variable_class.def("__add__", &Python::add<PV,PV,double>);
  polynomial_variable_class.def("__add__", &Python::add<PV,PV,R>);
  polynomial_variable_class.def("__add__", &Python::add<PV,PV,A>);
  polynomial_variable_class.def("__radd__", &Python::radd<PV,PV,double>);
  polynomial_variable_class.def("__radd__", &Python::radd<PV,PV,R>);
  polynomial_variable_class.def("__radd__", &Python::radd<PV,PV,A>);
  polynomial_variable_class.def("__sub__", &Python::sub<PV,PV,PV>);
  polynomial_variable_class.def("__sub__", &Python::sub<PV,PV,double>);
  polynomial_variable_class.def("__sub__", &Python::sub<PV,PV,R>);
  polynomial_variable_class.def("__sub__", &Python::sub<PV,PV,A>);
  polynomial_variable_class.def("__rsub__", &Python::rsub<PV,PV,double>);
  polynomial_variable_class.def("__rsub__", &Python::rsub<PV,PV,R>);
  polynomial_variable_class.def("__rsub__", &Python::rsub<PV,PV,A>);
  polynomial_variable_class.def("__mul__", &Python::mul<PV,PV,PV>);
  polynomial_variable_class.def("__mul__", &Python::mul<PV,PV,double>);
  polynomial_variable_class.def("__mul__", &Python::mul<PV,PV,R>);
  polynomial_variable_class.def("__mul__", &Python::mul<PV,PV,A>);
  polynomial_variable_class.def("__rmul__", &Python::rmul<PV,PV,double>);
  polynomial_variable_class.def("__rmul__", &Python::rmul<PV,PV,R>);
  polynomial_variable_class.def("__rmul__", &Python::rmul<PV,PV,A>);
  polynomial_variable_class.def("__div__", &Python::div<PV,PV,PV>);
  polynomial_variable_class.def("__div__", &Python::div<PV,PV,double>);
  polynomial_variable_class.def("__div__", &Python::div<PV,PV,R>);
  polynomial_variable_class.def("__div__", &Python::div<PV,PV,A>);
  polynomial_variable_class.def("__rdiv__", &Python::rdiv<PV,PV,double>);
  polynomial_variable_class.def("__rdiv__", &Python::rdiv<PV,PV,R>);
  polynomial_variable_class.def("__rdiv__", &Python::rdiv<PV,PV,A>);
  polynomial_variable_class.def("__pow__", &Python::pow<PV,PV,int>);
  polynomial_variable_class.def(self_ns::str(self));
  
  def("polynomial_constant",(PV(*)(size_type, smoothness_type, const double&))&PV::constant);
  def("polynomial_variable",(PV(*)(size_type, smoothness_type, const double&, size_type))&PV::variable);

  def("polynomial_constant",(PV(*)(size_type, smoothness_type, const A&))&PV::constant);
  def("polynomial_variable",(PV(*)(size_type, smoothness_type, const A&, size_type))&PV::variable);

  def("compose",(PV(*)(const PV&,const PV&))&compose);

  def("max",(PV(*)(const PV&,const PV&))&max);
  def("min",(PV(*)(const PV&,const PV&))&min);
  def("abs",(PV(*)(const PV&))&abs);
  def("neg",(PV(*)(const PV&))&neg);
  def("rec",(PV(*)(const PV&))&rec);
  def("pow",(PV(*)(const PV&, int))&pow);

  def("sqrt", (PV(*)(const PV&))&sqrt);
  def("exp", (PV(*)(const PV&))&exp);
  def("log", (PV(*)(const PV&))&log);
  def("sin", (PV(*)(const PV&))&sin);
  def("cos", (PV(*)(const PV&))&cos);
  def("tan", (PV(*)(const PV&))&tan);
  def("asin", (PV(*)(const PV&))&asin);
  def("acos", (PV(*)(const PV&))&acos);
  def("atan", (PV(*)(const PV&))&atan);

  def("rec_polynomial", (PV(*)(smoothness_type, const A&))&PV::rec);
  def("pow_polynomial",(PV(*)(smoothness_type, const A&, const uint&))&PV::pow);
  def("sqrt_polynomial", (PV(*)(smoothness_type, const A&))&PV::sqrt);
  def("exp_polynomial", (PV(*)(smoothness_type, const A&))&PV::exp);
  def("log_polynomial", (PV(*)(smoothness_type, const A&))&PV::log);
  def("sin_polynomial", (PV(*)(smoothness_type, const A&))&PV::sin);
  def("cos_polynomial", (PV(*)(smoothness_type, const A&))&PV::cos);
  def("tan_polynomial", (PV(*)(smoothness_type, const A&))&PV::tan);
  def("asin_polynomial", (PV(*)(smoothness_type, const A&))&PV::asin);
  def("acos_polynomial", (PV(*)(smoothness_type, const A&))&PV::acos);
  def("atan_polynomial", (PV(*)(smoothness_type, const A&))&PV::atan);

}

template void export_polynomial_variable<FloatPy>();
