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
#include "function/taylor_model.h"
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
  typedef typename traits<R>::arithmetic_type A;
  typedef typename traits<R>::interval_type I;

  class_< TaylorModel<R> > taylor_model_class("TaylorModel",init< Vector<I>, Vector<R>, const FunctionInterface<R>&, smoothness_type, smoothness_type>());
  taylor_model_class.def(init< size_type, size_type, smoothness_type, smoothness_type >());
  taylor_model_class.def(init< Vector<I>, Vector<R>, TaylorDerivative<I>, TaylorDerivative<I> >());
  taylor_model_class.def(init< TaylorModel<R> >());
  taylor_model_class.def("__call__",(Vector<I>(TaylorModel<R>::*)(const Vector<R>&)const) &TaylorModel<R>::evaluate);
  taylor_model_class.def("__call__",(Vector<I>(TaylorModel<R>::*)(const Vector<I>&)const) &TaylorModel<R>::evaluate);
  taylor_model_class.def("result_size", &TaylorModel<R>::result_size);
  taylor_model_class.def("argument_size", &TaylorModel<R>::argument_size);
  taylor_model_class.def("order", &TaylorModel<R>::order);
  taylor_model_class.def("smoothness", &TaylorModel<R>::smoothness);
  taylor_model_class.def("domain", &TaylorModel<R>::domain);
  taylor_model_class.def("centre", &TaylorModel<R>::centre);
  taylor_model_class.def("range", &TaylorModel<R>::range);
  //taylor_model_class.def("set",(void(TaylorModel<R>::*)(size_type,const MultiIndex&,const R&)) &TaylorModel<R>::set);
  //taylor_model_class.def("get",(R(TaylorModel<R>::*)(size_type,const MultiIndex&)const) &TaylorModel<R>::get);
  taylor_model_class.def("truncate",&TaylorModel<R>::truncate);
  taylor_model_class.def("evaluate",(Vector<I>(TaylorModel<R>::*)(const Vector<I>&)const) &TaylorModel<R>::evaluate);
  taylor_model_class.def("jacobian",(Matrix<I>(TaylorModel<R>::*)(const Vector<I>&)const) &TaylorModel<R>::jacobian);
  taylor_model_class.def("__add__",(TaylorModel<R>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &Ariadne::operator+<R>);
  taylor_model_class.def("__sub__",(TaylorModel<R>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &Ariadne::operator-<R>);
  taylor_model_class.def(self_ns::str(self));
 
  def("evaluate",(Vector<I>(TaylorModel<R>::*)(const Vector<R>&)const) &TaylorModel<R>::evaluate);
  def("evaluate",(Vector<I>(TaylorModel<R>::*)(const Vector<I>&)const) &TaylorModel<R>::evaluate);
  def("compose",(TaylorModel<R>(*)(const TaylorModel<R>&,const TaylorModel<R>&)) &compose);
  def("inverse",(TaylorModel<R>(*)(const TaylorModel<R>&)) &inverse);
  def("implicit",(TaylorModel<R>(*)(const TaylorModel<R>&)) &implicit);
  def("derivative",(TaylorModel<R>(*)(const TaylorModel<R>&, size_type)) &derivative);
}

template void export_taylor_model<FloatPy>();
