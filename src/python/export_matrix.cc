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

#include "base/numerical_type.h"
#include "base/interval.h"

#include "linear_algebra/matrix.h"
#include "linear_algebra/interval_matrix.h"

#include <boost/python.hpp>

#include "python/typedefs.h"
#include "python/python_utilities.h"

#include <utility>  //for std::pair

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

using boost::python::tuple;
using boost::python::extract;

inline Rational qmatrix_getitem(const QMatrix& A, tuple index) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  return A(i,j);
}

inline void qmatrix_setitem(QMatrix& A, tuple index, Rational x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=x;
}

inline void qmatrix_setitem_from_real(QMatrix& A, tuple index, Real x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=Ariadne::convert_to<Rational>(x);
}

inline void qmatrix_setitem_from_double(QMatrix& A, tuple index, double x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=Ariadne::convert_to<Rational>(x);
}

inline Real rmatrix_getitem(const RMatrix& A, tuple index) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  return A(i,j);
}

inline void rmatrix_setitem(RMatrix& A, tuple index, Real x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=x;
}

inline void rmatrix_setitem_from_double(RMatrix& A, tuple index, double x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=Ariadne::convert_to<Real>(x);
}

inline RVector rmatrix_vprod(const RMatrix& A, const RVector v) {
  return prod(A,v);
}

inline RInterval imatrix_getitem(const RIntervalMatrix& A, tuple index) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  return A(i,j);
}

inline void imatrix_setitem(RIntervalMatrix& A, tuple index, RInterval x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=x;
}

inline FMatrix rmatrix_inverse(const RMatrix& A) {
  return Ariadne::LinearAlgebra::inverse(A);
}

void export_matrix() {
  class_<RMatrix>("Matrix",init<int,int>())
    .def(init<RMatrix>())
    .def("__mul__",&rmatrix_vprod)
    .def("__getitem__",&rmatrix_getitem)
    .def("__setitem__",&rmatrix_setitem)
    .def("__setitem__",&rmatrix_setitem_from_double)
    .def(str(self))    // __str__
  ;
    
  def("inverse",&rmatrix_inverse);
  def("exp_approx",&Ariadne::LinearAlgebra::exp_approx<Real>);

  class_<QMatrix>("QMatrix",init<int,int>())
    .def(init<QMatrix>())
    .def("__getitem__",&qmatrix_getitem)
    .def("__setitem__",&qmatrix_setitem)
    .def("__setitem__",&qmatrix_setitem_from_real)
    .def("__setitem__",&qmatrix_setitem_from_double)
    .def(str(self))    // __str__
  ;
}

void export_interval_matrix() {
  class_<RIntervalMatrix>("IntervalMatrix",init<int,int>())
    .def(init<RIntervalMatrix>())
    .def("__getitem__",&imatrix_getitem)
    .def("__setitem__",&imatrix_setitem)
    .def(str(self))    // __str__
  ;
}
