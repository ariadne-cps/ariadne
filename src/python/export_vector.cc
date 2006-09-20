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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include "numeric/numerical_types.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;

#include <boost/python.hpp>
using namespace boost::python;

template<class R> 
inline
R
vector_get_item(const Vector<R>& v, int n) {
  if(n<0) {
    n+=v.size();
  }
  assert(0<=n);
  size_t m=size_t(n);
  assert(m<v.size());
  return v(m);
}

template<class R, class A> 
inline
void
vector_set_item(Vector<R>& v, int n, const A& x) {
  if(n<0) {
    n+=v.size();
  }
  assert(0<=n);
  size_t m=size_t(n);
  assert(m<v.size());
  R& r=v(m);
  r=R(x);
}


template<class R>
void export_vector()
{
  typedef Vector<R> Vec;
  
  class_<Vec>(python_name<R>("Vector").c_str(),init<uint>())
    .def(init<std::string>())
    .def(init<Vec>())
    .def("__len__", &Vec::size)
    .def("__getitem__",&vector_get_item<R>)
    .def("__setitem__",&vector_set_item<R,R>)
    .def("__setitem__",&vector_set_item<R,double>)
    .def("__add__",&add<Vec,Vec,Vec>)
    .def("__sub__",&sub<Vec,Vec,Vec>)
    .def("__rmul__",&mul<Vec,R,Vec>)
    .def("__mul__",&mul<Vec,Vec,R>)
    .def("__div__",&div<Vec,Vec,R>)
    .def(self_ns::str(self))
  ;
}


template<class R>
void export_interval_vector() {
  typedef Interval<R> Ivl;
  typedef Vector<R> Vec;
  typedef Vector< Interval<R> > IVec;
  
  class_<IVec>(python_name<R>("IntervalVector").c_str(),init<uint>())
    .def(init<std::string>())
    .def(init<IVec>())
    .def("__len__", &IVec::size)
    .def("__getitem__",&vector_get_item<Ivl>) 
    .def("__setitem__",&vector_set_item<Ivl,Ivl>)
    .def("__setitem__",&vector_set_item<Ivl,R>)
    .def("__setitem__",&vector_set_item<Ivl,double>)
    .def("__add__",&add<IVec,IVec,IVec>)
    .def("__add__",&add<IVec,IVec,Vec>)
    .def("__radd__",&add<IVec,Vec,IVec>)
    .def("__sub__",&sub<IVec,IVec,IVec>)
    .def("__sub__",&sub<IVec,IVec,Vec>)
    .def("__rsub__",&sub<IVec,Vec,IVec>)
    .def("__mul__",&mul<IVec,IVec,Ivl>)
    .def("__mul__",&mul<IVec,IVec,R>)
    .def("__rmul__",&mul<IVec,Ivl,IVec>)
    .def("__rmul__",&mul<IVec,R,IVec>)
    .def("__div__",&mul<IVec,IVec,Ivl>)
    .def("__div__",&mul<IVec,IVec,R>)
   .def(self_ns::str(self))
  ;
}

template void export_vector<Float64>();
template void export_vector<MPFloat>();
template void export_vector<Rational>();

template void export_interval_vector<Float64>();
template void export_interval_vector<MPFloat>();
