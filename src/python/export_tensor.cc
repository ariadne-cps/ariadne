/***************************************************************************
 *            python/export_tensor.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include <utility>  //for std::pair


#include "base/array.h"
#include "numeric/rational.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/tensor.h"

#include "python/python_utilities.h"
#include "python/python_float.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class R> 
inline R tensor_get_item(const Tensor<R>& T, boost::python::object obj) {
  boost::python::tuple index=extract<tuple>(obj);
  SizeArray ia(len(index));
  for(size_type i=0; i!=ia.size(); ++i) {
    ia[i]=extract<uint>(index[i]);
  }
  return T(ia);
}

template<class R, class A> 
inline void tensor_set_item(Tensor<R>& T, boost::python::object obj, const A& x) {
  boost::python::tuple index=extract<tuple>(obj);
  SizeArray ia(len(index));
  for(size_type i=0; i!=ia.size(); ++i) {
    ia[i]=extract<uint>(index[i]);
  }
  T(ia)=R(x);
}



template<class R>
void export_tensor()
{
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Tensor<R> Tns;
  
  class_<Tns>(python_name<R>("Tensor").c_str(),init<SizeArray>())
    //.def(init<std::string>())
    .def(init<Tns>())
    .def("__getitem__",&tensor_get_item<R>)
    .def("__setitem__",&tensor_set_item<R,R>)
    .def("__setitem__",&tensor_set_item<R,double>)
    .def(self_ns::str(self)) 
  ;
}

template<>
void export_tensor<Rational>() 
{
  typedef Rational R;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Tensor<R> Tns;
  
  class_<Tns>(python_name<R>("Tensor").c_str(),init<SizeArray>())
    //.def(init<std::string>())
    .def(init<Tns>())
    .def("__getitem__",&tensor_get_item<R>)
    .def("__setitem__",&tensor_set_item<R,R>)
    .def("__setitem__",&tensor_set_item<R,double>)
//    .def("__add__",&add<Tns,Tns,Tns>)
//    .def("__sub__",&sub<Tns,Tns,Tns>)
//    .def("__rmul__",&mul<Tns,R,Tns>)
//    .def("__mul__",&mul<Tns,Tns,R>)
//    .def("__div__",&mul<Tns,Tns,R>)
    .def(self_ns::str(self)) 
  ;
}

template<class R>
void export_interval_tensor() 
{
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Tensor<R> Tns;
  typedef Vector< Interval<R> > IVec;
  typedef Matrix< Interval<R> > IMx;
  typedef Tensor< Interval<R> > ITns;
  
  class_<ITns>(python_name<R>("FuzzyTensor").c_str(),init<SizeArray>())
    //.def(init<std::string>())
    .def(init<Tns>())
    .def(init<ITns>())
    .def("__getitem__",&tensor_get_item<I>)
    .def("__setitem__",&tensor_set_item<I,I>)
    .def("__setitem__",&tensor_set_item<I,R>)
    .def("__setitem__",&tensor_set_item<I,double>)
    .def(self_ns::str(self))
  ;
}

template void export_tensor<Float>();
template void export_tensor<Rational>();

template void export_interval_tensor<Float>();
