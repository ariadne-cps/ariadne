/***************************************************************************
 *            python/export_vector.cc
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
#include "linear_algebra/vector.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"

typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::LinearAlgebra::vector<Real> RVector;
typedef Ariadne::LinearAlgebra::vector< Ariadne::Interval<Real> > IVector;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

inline Real rvector_getitem(const RVector& v, uint i) {
  return v(i);
}

inline void rvector_setitem(RVector& v, uint i, Real x) {
  v(i)=x;
}

inline void rvector_setitem_from_double(RVector& v, uint i, double x) {
  v(i)=Ariadne::convert_to<Real>(x);
}

inline RInterval ivector_getitem(const IVector& v, uint i) {
  return v(i);
}

inline void ivector_setitem(IVector& v, uint i, RInterval x) {
  v(i)=x;
}


void export_vector() {
  class_<RVector>("Vector",init<int>())
    .def(init<RVector>())
    .def("__len__", &RVector::size)
    .def("__getitem__",&rvector_getitem)
    .def("__setitem__",&rvector_setitem)
    .def("__setitem__",&rvector_setitem_from_double)
    .def(str(self))    // __str__
  ;
}

void export_interval_vector() {
  class_<IVector>("IntervalVector",init<int>())
    .def(init<IVector>())
    .def("__len__", &IVector::size)
    .def("__getitem__",&ivector_getitem)
    .def("__setitem__",&ivector_setitem)
    .def(str(self))    // __str__
  ;
}
