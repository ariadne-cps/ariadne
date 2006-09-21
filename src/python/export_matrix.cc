/***************************************************************************
 *            python/export_matrix.cc
 *
 *  17 November 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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


#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;

#include <boost/python.hpp>
using namespace boost::python;

template<typename R> 
inline R matrix_get_item(const Matrix<R>& M, tuple index) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  return M(i,j);
}

template<typename R, typename A> 
inline void matrix_set_item(Matrix<R>& M, tuple index, const A& x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  M(i,j)=R(x);
}


template<typename R>
inline Matrix<R> matrix_inverse(const Matrix<R>& A) {
  return LinearAlgebra::inverse(A);
}

template<typename R>
void export_matrix()
{
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  
  class_<Mx>(python_name<R>("Matrix").c_str(),init<int,int>())
    .def(init<std::string>())
    .def(init<Mx>())
    .def("__getitem__",&matrix_get_item<R>)
    .def("__setitem__",&matrix_set_item<R,R>)
    .def("__setitem__",&matrix_set_item<R,double>)
    .def(self_ns::str(self)) 
  ;
}

template<>
void export_matrix<Rational>() 
{
  typedef Rational R;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  
  class_<Mx>(python_name<R>("Matrix").c_str(),init<int,int>())
    .def(init<std::string>())
    .def(init<Mx>())
    .def("__getitem__",&matrix_get_item<R>)
    .def("__setitem__",&matrix_set_item<R,R>)
    .def("__setitem__",&matrix_set_item<R,double>)
    .def("__add__",&add<Mx,Mx,Mx>)
    .def("__sub__",&sub<Mx,Mx,Mx>)
    .def("__rmul__",&mul<Mx,R,Mx>)
    .def("__mul__",&mul<Mx,Mx,R>)
    .def("__mul__",&mul<Vec,Mx,Vec>)
    .def("__mul__",&mul<Mx,Mx,Mx>)
    .def("__div__",&mul<Mx,Mx,R>)
    .def(self_ns::str(self)) 
  ;
}

template<typename R>
void export_interval_matrix() 
{
  typedef Interval<R> Ivl;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Vector< Interval<R> > IVec;
  typedef Matrix< Interval<R> > IMx;
  
  class_<IMx>(python_name<R>("IntervalMatrix").c_str(),init<int,int>())
    .def(init<std::string>())
    .def(init<IMx>())
    .def("__getitem__",&matrix_get_item<Ivl>)
    .def("__setitem__",&matrix_set_item<Ivl,Ivl>)
    .def("__setitem__",&matrix_set_item<Ivl,R>)
    .def("__setitem__",&matrix_set_item<Ivl,double>)
    .def(self_ns::str(self))
  ;
  
  def("inverse",&matrix_inverse<Ivl>);
  def("exp",&LinearAlgebra::exp<R>);
}

template void export_matrix<Float64>();
template void export_matrix<MPFloat>();
template void export_matrix<Rational>();

template void export_interval_matrix<Float64>();
template void export_interval_matrix<MPFloat>();
