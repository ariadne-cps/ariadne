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
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class R>  
Matrix<R> 
extract_matrix(boost::python::object obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  boost::python::list elements=extract<list>(obj);
  int m=boost::python::len(elements);
  list row=extract<list>(elements[0]);
  int n=boost::python::len(row);
  Matrix<R> A(m,n);
  for(int i=0; i!=m; ++i) {
    row=extract<list>(elements[i]);
    if(boost::python::len(row)!=n) {
      throw std::runtime_error("Matrix with rows of different sizes");
    }
    for(int j=0; j!=n; ++j) {
      extract<double> x(row[j]);
      if (x.check()) {
        A(i,j)=static_cast<R>(x);
      } else {
        extract<R> x(row[j]);
        A(i,j)=x;
      }
    }
  }
  return A;
}

template<class R> inline 
R matrix_get_item(const Matrix<R>& M, boost::python::object obj) {
  tuple index=extract<tuple>(obj);
  int i=extract<int>(index[0]);
  int j=extract<int>(index[1]);
  return M(i,j);
}

template<class R, class A> inline 
void matrix_set_item(Matrix<R>& M, boost::python::object obj, const A& x) {
  tuple index=extract<tuple>(obj);
  int i=extract<int>(index[0]);
  int j=extract<int>(index[1]);
  M(i,j)=R(x);
}


template<class R1,class R2> inline 
Vector<typename Numeric::traits<R1,R2>::arithmetic_type> 
matrix_vector_solve(const Matrix<R1>& A, const Vector<R2>& v) {
 return Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>(A).solve(v);
}

template<class R> inline 
Matrix<typename Numeric::traits<R>::arithmetic_type> 
matrix_inverse(const Matrix<R>& A) {
  return LinearAlgebra::inverse(A);
}

template<class R> inline 
Matrix<R> 
matrix_transpose(const Matrix<R>& A) {
  return A.transpose();
}

template<class R>
void export_matrix()
{
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Vector< Interval<R> > IVec;
  typedef Matrix<R> Mx;
  typedef Matrix< Interval<R> > IMx;

  class_<Mx>(python_name<R>("Matrix").c_str(),init<int,int>())
    .def(init<std::string>())
    .def(init<Mx>())
    .def("__getitem__",&matrix_get_item<R>)
    .def("__setitem__",&matrix_set_item<R,R>)
    .def("__setitem__",&matrix_set_item<R,double>)
    .def("__neg__",&neg<Mx,Mx>)
    .def("__add__",&add<IMx,Mx,Mx>)
    .def("__add__",&add<IMx,Mx,IMx>)
    .def("__sub__",&sub<IMx,Mx,Mx>)
    .def("__sub__",&sub<IMx,Mx,IMx>)
    .def("__rmul__",&mul<IMx,Mx,int,Mx,R>)
    .def("__rmul__",&mul<IMx,Mx,double,Mx,R>)
    .def("__rmul__",&mul<IMx,Mx,R>)
    .def("__rmul__",&mul<IMx,Mx,I>)
    .def("__mul__",&mul<IMx,Mx,int,Mx,R>)
    .def("__mul__",&mul<IMx,Mx,double,Mx,R>)
    .def("__mul__",&mul<IMx,Mx,R>)
    .def("__mul__",&mul<IMx,Mx,I>)
    .def("__mul__",&mul<IVec,Mx,Vec>)
    .def("__mul__",&mul<IVec,Mx,IVec>)
    .def("__mul__",&mul<IMx,Mx,Mx>)
    .def("__mul__",&mul<IMx,Mx,IMx>)
    .def("__div__",&div<IMx,Mx,int,Mx,R>)
    .def("__div__",&div<IMx,Mx,double,Mx,R>)
    .def("__div__",&div<IMx,Mx,R>)
    .def("__div__",&div<IMx,Mx,I>)
    .def(self_ns::str(self)) 
  ;
  def("solve",&matrix_vector_solve<R,R>);
  def("solve",&matrix_vector_solve<R,I>);
  def("transpose",&matrix_transpose<R>);
  def("inverse",&matrix_inverse<R>);

  def("extract_matrix",&extract_matrix<R>,"Extract an Ariadne matrix from a Python list of lists");

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
    /*
    .def("__getitem__",&matrix_get_item<R>)
    .def("__setitem__",&matrix_set_item<R,R>)
    .def("__setitem__",&matrix_set_item<R,double>)
    */
    .def("__neg__",&neg<Mx,Mx>)
    .def("__add__",&add<Mx,Mx,Mx>)
    .def("__sub__",&sub<Mx,Mx,Mx>)
    .def("__rmul__",&mul<Mx,Mx,int,Mx,R>)
    .def("__rmul__",&mul<Mx,Mx,double,Mx,R>)
    .def("__rmul__",&mul<Mx,Mx,R>)
    .def("__mul__",&mul<Mx,Mx,int,Mx,R>)
    .def("__mul__",&mul<Mx,Mx,double,Mx,R>)
    .def("__mul__",&mul<Mx,Mx,R>)
    .def("__mul__",&mul<Vec,Mx,Vec>)
    .def("__mul__",&mul<Mx,Mx,Mx>)
    .def("__div__",&div<Mx,Mx,int,Mx,R>)
    .def("__div__",&div<Mx,Mx,double,Mx,R>)
    .def("__div__",&div<Mx,Mx,R>)
    .def(self_ns::str(self)) 
  ;
  
  def("solve",&matrix_vector_solve<R,R>);
  def("transpose",&matrix_transpose<R>);
  def("inverse",&matrix_inverse<R>);
}

template<class R>
void export_interval_matrix() 
{
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Vector< Interval<R> > IVec;
  typedef Matrix< Interval<R> > IMx;
  
  class_<IMx>(python_name<R>("IntervalMatrix").c_str(),init<int,int>())
    .def(init<std::string>())
    .def(init<Mx>())
    .def(init<IMx>())
    /*
    .def("__getitem__",&matrix_get_item<I>)
    .def("__setitem__",&matrix_set_item<I,I>)
    .def("__setitem__",&matrix_set_item<I,R>)
    .def("__setitem__",&matrix_set_item<I,double>)
    */
    .def("__neg__",&neg<IMx,IMx>)
    .def("__add__",&add<IMx,IMx,Mx>)
    .def("__add__",&add<IMx,IMx,IMx>)
    .def("__sub__",&sub<IMx,IMx,Mx>)
    .def("__sub__",&sub<IMx,IMx,IMx>)
    .def("__rmul__",&mul<IMx,IMx,int,IMx,R>)
    .def("__rmul__",&mul<IMx,IMx,double,IMx,R>)
    .def("__rmul__",&mul<IMx,IMx,R>)
    .def("__rmul__",&mul<IMx,IMx,I>)
    .def("__mul__",&mul<IMx,IMx,int,IMx,R>)
    .def("__mul__",&mul<IMx,IMx,double,IMx,R>)
    .def("__mul__",&mul<IMx,IMx,R>)
    .def("__mul__",&mul<IMx,IMx,I>)
    .def("__mul__",&mul<IVec,IMx,Vec>)
    .def("__mul__",&mul<IVec,IMx,IVec>)
    .def("__mul__",&mul<IMx,IMx,Mx>)
    .def("__mul__",&mul<IMx,IMx,IMx>)
    .def("__div__",&div<IMx,IMx,int,IMx,R>)
    .def("__div__",&div<IMx,IMx,double,IMx,R>)
    .def("__div__",&div<IMx,IMx,R>)
    .def("__div__",&div<IMx,IMx,I>)
    .def(self_ns::str(self))
  ;
  
  def("solve",&matrix_vector_solve<R,I>);
  def("solve",&matrix_vector_solve<I,R>);
  def("solve",&matrix_vector_solve<I,I>);
  def("transpose",&matrix_transpose<I>);
  def("inverse",&matrix_inverse<I>);
  def("exp",&LinearAlgebra::exp<R>);
}

template void export_matrix<Float64>();
template void export_matrix<MPFloat>();
template void export_matrix<Rational>();

template void export_interval_matrix<Float64>();
template void export_interval_matrix<MPFloat>();
