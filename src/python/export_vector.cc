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
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class R>  
Vector<R> 
extract_vector(boost::python::object obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  list elements=extract<list>(obj);
  int n=boost::python::len(elements);
  Vector<R> v(n);
  for(int i=0; i!=n; ++i) {
    extract<double> x(elements[i]);
    if (x.check()) {
      v(i)=static_cast<R>(x);
    } else {
      extract<R> x(elements[i]);
      v(i)=x;
    }
  }
  return v;
}


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
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Vector< Interval<R> > IVec;
  
  class_<Vec>(python_name<R>("Vector").c_str(),init<uint>())
    .def(init<std::string>())
    .def(init<Vec>())
    .def("__len__", &Vec::size)
    .def("__getitem__",&vector_get_item<R>)
    .def("__setitem__",&vector_set_item<R,R>)
    .def("__setitem__",&vector_set_item<R,double>)
    .def("__neg__",&neg<Vec,Vec>)
    .def("__add__",&add<IVec,Vec,Vec>)
    .def("__add__",&add<IVec,Vec,IVec>)
    .def("__sub__",&sub<IVec,Vec,Vec>)
    .def("__sub__",&sub<IVec,Vec,IVec>)
    .def("__rmul__",&mul<IVec,Vec,int,Vec,R>)
    .def("__rmul__",&mul<IVec,Vec,double,Vec,R>)
    .def("__rmul__",&mul<IVec,Vec,R>)
    .def("__rmul__",&mul<IVec,Vec,I>)
    .def("__mul__",&mul<IVec,Vec,int,Vec,R>)
    .def("__mul__",&mul<IVec,Vec,double,Vec,R>)
    .def("__mul__",&mul<IVec,Vec,R>)
    .def("__mul__",&mul<IVec,Vec,I>)
    .def("__div__",&div<IVec,Vec,int,Vec,R>)
    .def("__div__",&div<IVec,Vec,double,Vec,R>)
    .def("__div__",&div<IVec,Vec,R>)
    .def("__div__",&div<IVec,Vec,I>)
    .def(self_ns::str(self))
  ;

  def("extract_vector",&extract_vector<R>,"Extract an Ariadne vector from a Python list");
}

template<>
void export_vector<Rational>()
{
  typedef Rational R;
  typedef Vector<R> Vec;
  
  class_<Vec>(python_name<R>("Vector").c_str(),init<uint>())
    .def(init<std::string>())
    .def(init<Vec>())
    .def("__len__", &Vec::size)
    .def("__getitem__",&vector_get_item<R>)
    .def("__setitem__",&vector_set_item<R,R>)
    .def("__setitem__",&vector_set_item<R,double>)
    .def("__neg__",&neg<Vec,Vec>)
    .def("__add__",&add<Vec,Vec,Vec>)
    .def("__sub__",&sub<Vec,Vec,Vec>)
    .def("__rmul__",&mul<Vec,Vec,int,Vec,R>)
    .def("__rmul__",&mul<Vec,Vec,double,Vec,R>)
    .def("__rmul__",&mul<Vec,Vec,R,Vec,R>)
    .def("__mul__",&mul<Vec,Vec,int,Vec,R>)
    .def("__mul__",&mul<Vec,Vec,double,Vec,R>)
    .def("__mul__",&mul<Vec,Vec,R,Vec,R>)
    .def("__div__",&div<Vec,Vec,int,Vec,R>)
    .def("__div__",&div<Vec,Vec,double,Vec,R>)
    .def("__div__",&div<Vec,Vec,R,Vec,R>)
    .def(self_ns::str(self))
  ;
}


template<class R>
void export_interval_vector() {
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Vector< Interval<R> > IVec;
  
  class_<IVec>(python_name<R>("IntervalVector").c_str(),init<uint>())
    .def(init<std::string>())
    .def(init<Vec>())
    .def(init<IVec>())
    .def("__len__", &IVec::size)
    .def("__getitem__",&vector_get_item<I>) 
    .def("__setitem__",&vector_set_item<I,I>)
    .def("__setitem__",&vector_set_item<I,R>)
    .def("__setitem__",&vector_set_item<I,double>)
    .def("__neg__",&neg<IVec,IVec>)
    .def("__add__",&add<IVec,IVec,Vec>)
    .def("__add__",&add<IVec,IVec,IVec>)
    .def("__sub__",&sub<IVec,IVec,Vec>)
    .def("__sub__",&sub<IVec,IVec,IVec>)
    .def("__rmul__",&mul<IVec,IVec,int,IVec,R>)
    .def("__rmul__",&mul<IVec,IVec,double,IVec,R>)
    .def("__rmul__",&mul<IVec,IVec,R>)
    .def("__rmul__",&mul<IVec,IVec,I>)
    .def("__mul__",&mul<IVec,IVec,int,IVec,R>)
    .def("__mul__",&mul<IVec,IVec,double,IVec,R>)
    .def("__mul__",&mul<IVec,IVec,R>)
    .def("__mul__",&mul<IVec,IVec,I>)
    .def("__div__",&div<IVec,IVec,int,IVec,R>)
    .def("__div__",&div<IVec,IVec,double,IVec,R>)
    .def("__div__",&div<IVec,IVec,R>)
    .def("__div__",&div<IVec,IVec,I>)
    .def(self_ns::str(self))
  ;
}

template void export_vector<Float64>();
template void export_vector<MPFloat>();
template void export_vector<Rational>();

template void export_interval_vector<Float64>();
template void export_interval_vector<MPFloat>();
