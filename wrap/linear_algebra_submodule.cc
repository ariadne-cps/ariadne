/***************************************************************************
 *            linear_algebra_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 
#include "config.h"

#include "numeric.h"
#include "vector.h"
#include "matrix.h"

#include "utilities.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;

namespace Ariadne {

template<class X> const char* python_name(const char* name);
template<> const char* python_name<Float>(const char* name) { 
    return (std::string("")+name).c_str(); }
template<> const char* python_name<Interval>(const char* name) { 
    return (std::string("I")+name).c_str(); }
template<> const char* python_name<Rational>(const char* name) { 
    return (std::string("Q")+name).c_str(); }


template<class X> 
void read(Vector<X>& v, const boost::python::object& obj) 
{
  array<X> a; 
  read_array(a,obj);
  v=Vector<X>(a.size(),a.begin());
}

template<class X> 
void read(Matrix<X>& A, const boost::python::object& obj) 
{
  boost::python::list elements=boost::python::extract<boost::python::list>(obj);
  int m=boost::python::len(elements);
  boost::python::list row=boost::python::extract<boost::python::list>(elements[0]);
  int n=boost::python::len(row);
  A=Matrix<X>(m,n);
  for(int i=0; i!=m; ++i) {
    row=boost::python::extract<boost::python::list>(elements[i]);
    if(boost::python::len(row)!=n) { throw std::runtime_error("Matrix with rows of different sizes"); }
    for(int j=0; j!=n; ++j) { X x; read_scalar(x,row[j]); A.set(i,j,x); }
  }
}


template<class X> Vector<X>*
make_vector(const boost::python::object& obj) 
{
    Vector<X>* vptr=new Vector<X>();
    read(*vptr,obj);
    return vptr;
}

template<class X> Matrix<X>*
make_matrix(const boost::python::object& obj) 
{
    Matrix<X>* Aptr=new Matrix<X>();
    read(*Aptr,obj);
    return Aptr;
}

template<class X0, class X1>
Matrix<X0> __neg__(const Matrix<X1>& A1) {
    return -A1;
}

template<class X0, class X1, class X2>
Matrix<X0> __add__(const Matrix<X1>& A1, const X2& A2) {
    return A1+A2;
}

template<class X0, class X1, class X2>
Matrix<X0> __sub__(const Matrix<X1>& A1, const X2& A2) {
    return A1-A2;
}

template<class X0, class X1, class X2>
Matrix<X0> __mul__(const Matrix<X1>& A1, const X2& s2) {
    return A1*s2;
}

template<class X0, class X1, class X2>
Matrix<X0> __div__(const Matrix<X1>& A1, const X2& s2) {
    return A1/s2;
}

template<class X0, class X1, class X2>
Matrix<X0> __mul__(const X1& s1, const Matrix<X2>& A2) {
    return s1*A2;
}

template<class X0, class X1, class X2>
Vector<X0> __mvmul__(const Matrix<X1>& A1, const Vector<X2>& v2) {
    return prod(A1,v2); 
}

template<class X0, class X1, class X2>
Matrix<X0> __mmmul__(const Matrix<X1>& A1, const Matrix<X2>& A2) {
    return prod(A1,A2); 
}


}


template<class X>
void export_vector()
{
    typedef Vector<X> V;

    class_< Vector<X> > vector_class(python_name<X>("Vector"),init<int>());
    vector_class.def("__init__", make_constructor(&make_vector<X>) );
    vector_class.def("__len__", &Vector<X>::size);
    vector_class.def("__setitem__", (void(Vector<X>::*)(size_t,const double&)) &Vector<X>::set);
    vector_class.def("__setitem__", (void(Vector<X>::*)(size_t,const X&)) &Vector<X>::set);
    vector_class.def("__getitem__", &Vector<X>::get, return_value_policy<copy_const_reference>());
    vector_class.def(-self);
    vector_class.def(self + self);
    vector_class.def(self - self);
    vector_class.def(X() * self);
    vector_class.def(self * X());
    vector_class.def(self / X());
    vector_class.def(double() * self);
    vector_class.def(self * double());
    vector_class.def(self / double());
    vector_class.def(self += self);
    vector_class.def(self -= self);
    vector_class.def(self *= double());
    vector_class.def(self *= X());
    vector_class.def(self /= double());
    vector_class.def(self /= X());
    vector_class.def(boost::python::self_ns::str(self));

    def("norm",(X(*)(const Vector<X>&)) &norm);
}


template<class X>
void export_matrix()
{
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),init<size_t,size_t>());
    matrix_class.def("__init__", make_constructor(&make_matrix<X>) );
    matrix_class.def("rows", &Matrix<X>::row_size);
    matrix_class.def("columns", &Matrix<X>::column_size);
    matrix_class.def("row_size", &Matrix<X>::row_size);
    matrix_class.def("column_size", &Matrix<X>::column_size);
    matrix_class.def("__setitem__", (void(Matrix<X>::*)(size_t,size_t,const double&)) &Matrix<X>::set);
    matrix_class.def("__setitem__", (void(Matrix<X>::*)(size_t,size_t,const X&)) &Matrix<X>::set);
    matrix_class.def("__getitem__", &Matrix<X>::get, return_value_policy<copy_const_reference>());
    matrix_class.def("__neg__", (Matrix<X>(*)(const Matrix<X>&)) &__neg__);
    matrix_class.def("__add__", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &__add__);
    matrix_class.def("__sub__", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &__sub__);
    matrix_class.def("__mul__", (Matrix<X>(*)(const Matrix<X>&,const X&)) &__mul__);
    matrix_class.def("__rmul__", (Matrix<X>(*)(const X&,const Matrix<X>&)) &__mul__);
    matrix_class.def("__div__", (Matrix<X>(*)(const Matrix<X>&,const X&)) &__div__);
    matrix_class.def("__mul__",(Vector<X>(*)(const Matrix<X>&,const Vector<X>&))&__mvmul__);
    matrix_class.def("__mul__",(Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&))&__mmmul__);
    //matrix_class.def("inverse", &inverse<X>);
    //matrix_class.def("determinant", &Matrix::determinant);
    //matrix_class.def("transpose", &Matrix::transpose);
    //matrix_class.def("solve", &Matrix::solve);
    matrix_class.def(boost::python::self_ns::str(self));    // __str__

    def("norm",(X(*)(const Matrix<X>&)) &norm);
    def("inverse",(Matrix<X>(*)(const Matrix<X>&)) &inverse);
}

template void export_vector<Float>();
template void export_vector<Interval>();
template void export_matrix<Float>();
template void export_matrix<Interval>();


#ifdef HAVE_GMPXX_H
template void export_vector<Rational>();
template void export_matrix<Rational>();
#endif


void linear_algebra_submodule() {
    export_vector<Float>();
    export_vector<Interval>();
    export_matrix<Float>();
    export_matrix<Interval>();
#ifdef HAVE_GMPXX_H
    export_vector<Rational>();
    export_matrix<Rational>();
#endif
}
