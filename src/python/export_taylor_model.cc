/***************************************************************************
 *            python/export_taylor_model.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/taylor_derivative.h"
#include "differentiation/sparse_differential.h"
#include "function/taylor_model.h"
#include "function/approximate_taylor_model.h"
#include "function/function_interface.h"
#include "geometry/box.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
TaylorModel<R> 
add(const TaylorModel<R>& tm1, const TaylorModel<R>& tm2) {
  return Ariadne::add(tm1,tm2);
}

template<class R>
TaylorModel<R> 
sub(const TaylorModel<R>& tm1, const TaylorModel<R>& tm2) {
  return Ariadne::sub(tm1,tm2);
}


template<class R>
void export_taylor_model() 
{
  typedef TaylorModel<R> Model;
  typedef typename traits<R>::arithmetic_type A;
  typedef typename traits<R>::interval_type I;

  class_< Model > function_model_class("TaylorModel",init< Vector<I>, Vector<R>, const FunctionInterface<R>&, smoothness_type, smoothness_type>());
  function_model_class.def(init< size_type, size_type, smoothness_type, smoothness_type >());
  //function_model_class.def(init< Vector<I>, Vector<R>, TaylorDerivative<I>, TaylorDerivative<I> >());
  function_model_class.def(init< Model >());
  //function_model_class.def("__call__",(Vector<I>(Model::*)(const Vector<R>&)const) &Model::evaluate);
  function_model_class.def("__call__",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  function_model_class.def("result_size", &Model::result_size);
  function_model_class.def("argument_size", &Model::argument_size);
  function_model_class.def("order", &Model::order);
  function_model_class.def("smoothness", &Model::smoothness);
  function_model_class.def("domain", &Model::domain);
  function_model_class.def("centre", &Model::centre);
  function_model_class.def("range", &Model::range);
  //function_model_class.def("set",(void(Model::*)(size_type,const MultiIndex&,const R&)) &Model::set);
  //function_model_class.def("get",(R(Model::*)(size_type,const MultiIndex&)const) &Model::get);
  function_model_class.def("truncate",&Model::truncate);
  function_model_class.def("evaluate",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  function_model_class.def("jacobian",(Matrix<I>(Model::*)(const Vector<I>&)const) &Model::jacobian);
  function_model_class.def("__add__",(Model(*)(const Model&,const Model&)) &Ariadne::operator+<R>);
  function_model_class.def("__sub__",(Model(*)(const Model&,const Model&)) &Ariadne::operator-<R>);
  function_model_class.def(self_ns::str(self));
 
  //def("evaluate",(Vector<I>(Model::*)(const Vector<R>&)const) &Model::evaluate);
  def("evaluate",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  def("compose",(Model(*)(const Model&,const Model&)) &compose);
  def("inverse",(Model(*)(const Model&)) &inverse);
  def("implicit",(Model(*)(const Model&)) &implicit);
  def("derivative",(Model(*)(const Model&, size_type)) &derivative);
}

template<class R>
void export_approximate_taylor_model() 
{
  typedef ApproximateTaylorModel<R> Model;
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef typename traits<R>::interval_type I;

  class_< Model > function_model_class("ApproximateTaylorModel",init< Vector<I>, Vector<R>, const FunctionInterface<R>&, smoothness_type, smoothness_type>());
  function_model_class.def(init< size_type, size_type, smoothness_type, smoothness_type >());
  function_model_class.def(init< Vector<I>, Vector<A>, SparseDifferentialVector<A> >());
  function_model_class.def(init< Vector<I>, Vector<R>, SparseDifferentialVector<A> >());
  function_model_class.def(init< Model >());
  //function_model_class.def("__call__",(Vector<I>(Model::*)(const Vector<R>&)const) &Model::evaluate);
  function_model_class.def("__call__",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  function_model_class.def("result_size", &Model::result_size);
  function_model_class.def("argument_size", &Model::argument_size);
  function_model_class.def("order", &Model::order);
  function_model_class.def("smoothness", &Model::smoothness);
  function_model_class.def("domain", &Model::domain);
  function_model_class.def("centre", &Model::centre);
  function_model_class.def("range", &Model::range);
  function_model_class.def("expansion", &Model::expansion, return_value_policy<copy_const_reference>());
  //function_model_class.def("set",(void(Model::*)(size_type,const MultiIndex&,const R&)) &Model::set);
  //function_model_class.def("get",(R(Model::*)(size_type,const MultiIndex&)const) &Model::get);
  function_model_class.def("truncate",&Model::truncate);
  function_model_class.def("evaluate",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  function_model_class.def("jacobian",(Matrix<I>(Model::*)(const Vector<I>&)const) &Model::jacobian);
  function_model_class.def("__add__",(Model(*)(const Model&,const Model&)) &Ariadne::operator+<R>);
  function_model_class.def("__sub__",(Model(*)(const Model&,const Model&)) &Ariadne::operator-<R>);
  function_model_class.def(self_ns::str(self));
 
  //def("evaluate",(Vector<I>(Model::*)(const Vector<R>&)const) &Model::evaluate);
  def("evaluate",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  def("project",(Model(*)(const Model&,const Slice&)) &project);
  def("join",(Model(*)(const Model&,const Model&)) &join);
  def("compose",(Model(*)(const Model&,const Model&)) &compose);
  def("inverse",(Model(*)(const Model&)) &inverse);
  def("implicit",(Model(*)(const Model&)) &implicit);
  def("derivative",(Model(*)(const Model&, size_type)) &derivative);
  def("flow",(Model(*)(const Model&)) &flow);
  //def("integrate",(Model(*)(const Model&,const R&)) &integrate);
  def("hitting",(Model(*)(const Model&,const Model&)) &hitting);
  def("solve",(Vector<I>(*)(const Model&,const Vector<R>&)) &solve);
}

template void export_taylor_model< FloatPy >();
template void export_approximate_taylor_model< FloatPy >();
