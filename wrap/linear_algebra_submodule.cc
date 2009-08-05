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
#include <boost/python/slice.hpp>

using namespace boost::python;

using namespace Ariadne;

namespace Ariadne {


#ifdef HAVE_GMPXX_H
template<> const char* python_name<Rational>(const char* name) {
    return (std::string("Q")+name).c_str(); }
#endif


template<class X>
X __vgetitem__(const Vector<X>& v, int i)
{
    if(i<0) { i+=v.size(); }
    ARIADNE_ASSERT(0<=i && uint(i)<v.size());
    return v[i];
}


template<class X>
Vector<X> __vgetslice__(const Vector<X>& v, int start, int stop)
{
    if(start<0) { start+=v.size(); }
    if(stop<0) { stop+=v.size(); }
    ARIADNE_ASSERT(0<=start && start<=stop && uint(stop)<=v.size());
    return project(v,range(start,stop));
}

/* This code is for start:stop:step slices, which we don't use
template<class X>
Vector<X> __vgetslice__(const Vector<X>& v, const boost::python::slice& slc)
{
    int start=boost::python::extract<int>(slc.start());
    int stop=boost::python::extract<int>(slc.stop());
    int step=boost::python::extract<int>(slc.step());
    ARIADNE_ASSERT(step==1);
    return project(v,range(start,stop));
}
*/

template<class X>
void read(Vector<X>& v, const boost::python::object& obj)
{
  array<X> a;
  read_list_array(a,obj);
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
    for(int j=0; j!=n; ++j) { X x; read(x,row[j]); A.set(i,j,x); }
  }
}

template void read(Matrix<double>&,const object&);

template<class X> Vector<X>*
make_vector(const boost::python::list& lst)
{
    Vector<X>* vptr=new Vector<X>(boost::python::len(lst));
    for(uint i=0; i!=vptr->size(); ++i) {
        (*vptr)[i]=boost::python::extract<X>(lst[i]);
    }
    return vptr;
}

template<class X> Matrix<X>*
make_matrix(const boost::python::list& lst)
{
    list const& elements=lst;
    int m=len(elements);
    list row=extract<list>(elements[0]);
    int n=len(row);
    Matrix<X>* Aptr=new Matrix<X>(m,n);
    for(int i=0; i!=m; ++i) {
        row=extract<list>(elements[i]);
        if(len(row)!=n) { throw std::runtime_error("Matrix with rows of different sizes"); }
        for(int j=0; j!=n; ++j) { Aptr->set(i,j,extract<X>(row[j])); }
    }
    //read(*Aptr,obj);
    return Aptr;
}


template<class X>
X __mgetitem__(const Matrix<X>& A, const boost::python::tuple& tup)
{
    uint i=boost::python::extract<uint>(tup[0]);
    uint j=boost::python::extract<uint>(tup[1]);
    return A[i][j];
}

template<class X, class XX>
void __msetitem__(Matrix<X>& A, const boost::python::tuple& tup, const XX& x)
{
    uint i=boost::python::extract<uint>(tup[0]);
    uint j=boost::python::extract<uint>(tup[1]);
    A[i][j]=X(x);
}




template<class X0, class X1>
Vector<X0> __vpos__(const Vector<X1>& v1) {
    return +v1;
}

template<class X0, class X1>
Vector<X0> __vneg__(const Vector<X1>& v1) {
    return -v1;
}

template<class X0, class X1, class X2>
Vector<X0> __vvadd__(const Vector<X1>& v1, const Vector<X2>& v2) {
    return v1+v2;
}

template<class X0, class X1, class X2>
Vector<X0> __vvsub__(const Vector<X1>& v1, const Vector<X2>& v2) {
    return v1-v2;
}

template<class X0, class X1, class X2>
Vector<X0> __svmul__(const X1& s1, const Vector<X2>& v2) {
    return s1*v2;
}

template<class X0, class X1, class X2>
Vector<X0> __vsmul__(const Vector<X1>& v1, const X2& s2) {
    return v1*s2;
}

template<class X0, class X1, class X2>
Vector<X0> __vsdiv__(const Vector<X1>& v1, const X2& s2) {
    return v1/s2;
}


template<class X0, class X1>
Matrix<X0> __mpos__(const Matrix<X1>& A1) {
    return +A1;
}

template<class X0, class X1>
Matrix<X0> __mneg__(const Matrix<X1>& A1) {
    return -A1;
}

template<class X0, class X1, class X2>
Matrix<X0> __mmadd__(const Matrix<X1>& A1, const Matrix<X2>& A2) {
    return A1+A2;
}

template<class X0, class X1, class X2>
Matrix<X0> __mmsub__(const Matrix<X1>& A1, const Matrix<X2>& A2) {
    return A1-A2;
}

template<class X0, class X1, class X2>
Matrix<X0> __smmul__(const X1& s1, const Matrix<X2>& A2) {
    return s1*A2;
}

template<class X0, class X1, class X2>
Matrix<X0> __msmul__(const Matrix<X1>& A1, const X2& s2) {
    return A1*s2;
}

template<class X0, class X1, class X2>
Matrix<X0> __msdiv__(const Matrix<X1>& A1, const X2& s2) {
    return A1/s2;
}

template<class X0, class X1, class X2>
Vector<X0> __mvmul__(const Matrix<X1>& A1, const Vector<X2>& v2) {
    return prod(A1,v2);
}

template<class X0, class X1, class X2>
Matrix<X0> __mmmul__(const Matrix<X1>& A1, const Matrix<X2>& A2) {
    return prod(A1,A2);
}


boost::python::tuple
wrap_orthogonal_decomposition(const Matrix<Float>& A)
{
    Ariadne::tuple<Matrix<Float>,Matrix<Float>, PivotMatrix > QRP=Ariadne::orthogonal_decomposition(A);
    return boost::python::make_tuple(QRP.first,QRP.second,pivot_matrix(QRP.third));
}

boost::python::tuple
wrap_triangular_factor(const Matrix<Float>& A)
{
    Ariadne::tuple<Matrix<Float>,PivotMatrix > RP=Ariadne::triangular_factor(A);
    return boost::python::make_tuple(RP.first,pivot_matrix(RP.second));
}

boost::python::tuple
wrap_triangular_decomposition(const Matrix<Float>& A)
{
    Ariadne::tuple<PivotMatrix,Matrix<Float>,Matrix<Float> > PLU=Ariadne::triangular_decomposition(A);
    return boost::python::make_tuple(pivot_matrix(PLU.first),PLU.second,PLU.third);
}




}


template<class X>
void export_vector_class(class_<Vector<X> >& vector_class)
{
    vector_class.def("__init__", make_constructor(&make_vector<X>) );
    vector_class.def(init<int>());
    vector_class.def(init<int,X>());
    vector_class.def("__len__", &Vector<X>::size);
    vector_class.def("__setitem__", (void(Vector<X>::*)(size_t,const double&)) &Vector<X>::set);
    vector_class.def("__setitem__", (void(Vector<X>::*)(size_t,const X&)) &Vector<X>::set);
    vector_class.def("__getitem__", &__vgetitem__<X>);
    vector_class.def("__getslice__", &__vgetslice__<X>);
    vector_class.def("__neg__", &__vneg__<X,X>);
    vector_class.def("__repr__",&__cstr__<Vector<X> >);
    vector_class.def("__pos__",__pos__< Vector<X>, Vector<X> >);
    vector_class.def("__neg__",__neg__< Vector<X>, Vector<X> >);
    vector_class.def(boost::python::self_ns::str(self));
    vector_class.def("unit",&Vector<X>::unit);
    vector_class.def("basis",&Vector<X>::basis);
    vector_class.staticmethod("unit");
    vector_class.staticmethod("basis");

    def("norm",(X(*)(const Vector<X>&)) &norm);
    def("join",(Vector<X>(*)(const Vector<X>&,const Vector<X>&)) &join);
    def("join",(Vector<X>(*)(const Vector<X>&,const X&)) &join);

    to_python_converter< array< Vector<X> >, array_to_python_list< Vector<X> > >();


}

template<class X, class Y>
void export_vector_conversion(class_<Vector<X> >& vector_class)
{
    vector_class.def(init<int,Y>());
    vector_class.def(init< Vector<Y> >());
}

template<class R, class X, class Y>
void export_vector_arithmetic(class_<Vector<X> >& vector_class)
{
    vector_class.def("__add__",__vvadd__<R,X,Y>);
    vector_class.def("__sub__",__vvsub__<R,X,Y>);
    vector_class.def("__rmul__",__vsmul__<R,X,Y>);
    vector_class.def("__mul__",__vsmul__<R,X,Y>);
    vector_class.def("__div__",__vsdiv__<R,X,Y>);
    //vector_class.def(self+=Vector<Y>());
    //vector_class.def(self-=Vector<Y>());
    //vector_class.def(self*=Y());
    //vector_class.def(self/=Y());
}


template<class X> void export_vector()
{
    class_< Vector<X> > vector_class(python_name<X>("Vector"),no_init);
    export_vector_class<X>(vector_class);
    export_vector_conversion<X,X>(vector_class);
    export_vector_conversion<X,Float>(vector_class);
    export_vector_arithmetic<X,X,X>(vector_class);
    export_vector_arithmetic<X,X,Float>(vector_class);
}


template<> void export_vector<Float>()
{
    class_< Vector<Float> > float_vector_class("Vector",no_init);
    export_vector_class<Float>(float_vector_class);
    export_vector_conversion<Float,Float>(float_vector_class);
    export_vector_arithmetic<Float,Float,Float>(float_vector_class);
    export_vector_arithmetic<Interval,Float,Interval>(float_vector_class);
}

template<> void export_vector<Interval>()
{
    class_< Vector<Interval> > interval_vector_class("IVector",no_init);
    export_vector_class<Interval>(interval_vector_class);
    export_vector_conversion<Interval,Float>(interval_vector_class);
    export_vector_conversion<Interval,Interval>(interval_vector_class);
    export_vector_arithmetic<Interval,Interval,Float>(interval_vector_class);
    export_vector_arithmetic<Interval,Interval,Interval>(interval_vector_class);
    def("midpoint", (Vector<Float>(*)(const Vector<Interval>&)) &midpoint);
    def("subset", (bool(*)(const Vector<Interval>&, const Vector<Interval>&)) &subset);
    implicitly_convertible< Vector<Float>, Vector<Interval> >();
}





template<class X>
void export_matrix_class(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def("__init__", make_constructor(&make_matrix<X>) );
    matrix_class.def(init<int,int>());
    matrix_class.def("rows", &Matrix<X>::row_size);
    matrix_class.def("columns", &Matrix<X>::column_size);
    matrix_class.def("row_size", &Matrix<X>::row_size);
    matrix_class.def("column_size", &Matrix<X>::column_size);
    matrix_class.def("__setitem__", &__msetitem__<X,double>);
    matrix_class.def("__setitem__", &__msetitem__<X,X>);
    matrix_class.def("__getitem__", &__mgetitem__<X>);
    matrix_class.def("__pos__", &__mpos__<X,X>);
    matrix_class.def("__neg__", &__mneg__<X,X>);
    //matrix_class.def("inverse", &inverse<X>);
    //matrix_class.def("determinant", &Matrix::determinant);
    //matrix_class.def("transpose", &Matrix::transpose);
    //matrix_class.def("solve", &Matrix::solve);
    matrix_class.def(boost::python::self_ns::str(self));    // __str__
    matrix_class.def("__repr__",&__cstr__<Matrix<X> >);

    def("norm",(X(*)(const Matrix<X>&)) &norm);
    def("transpose",(Matrix<X>(*)(const Matrix<X>&)) &transpose);

    def("inverse",(Matrix<X>(*)(const Matrix<X>&)) &inverse);
    def("solve",(Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &solve);
    def("solve",(Vector<X>(*)(const Matrix<X>&,const Vector<X>&)) &solve);

    def("triangular_decomposition",&wrap_triangular_decomposition);
    def("orthogonal_decomposition",&wrap_orthogonal_decomposition);
    def("triangular_factor",&wrap_triangular_factor);
    def("triangular_multiplier",(Matrix<X>(*)(const Matrix<X>&)) &triangular_multiplier);
    def("row_norms",(Vector<X>(*)(const Matrix<X>&)) &row_norms);
    def("normalise_rows",(Matrix<X>(*)(const Matrix<X>&)) &normalise_rows);
}

template<class X, class Y>
void export_matrix_conversion(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def(init< Matrix<Y> >());
}

template<class R, class X, class Y>
void export_matrix_arithmetic(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def("__add__", &__mmadd__<R,X,Y>);
    matrix_class.def("__sub__", &__mmsub__<R,X,Y>);
    matrix_class.def("__mul__", &__msmul__<R,X,Y>);
    matrix_class.def("__rmul__", &__msmul__<R,X,Y>);
    matrix_class.def("__div__", &__msdiv__<R,X,Y>);
    matrix_class.def("__mul__", &__mvmul__<R,X,Y>);
    matrix_class.def("__mul__", &__mmmul__<R,X,Y>);
}


template<class X> void export_matrix()
{
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),no_init);
    export_matrix_class<X>(matrix_class);
    export_matrix_conversion<X,Float>(matrix_class);
    export_matrix_conversion<X,X>(matrix_class);
    export_matrix_arithmetic<X,X,X>(matrix_class);
    export_matrix_arithmetic<X,X,Float>(matrix_class);
}

template<> void export_matrix<Float>()
{
    class_< Matrix<Float> > matrix_class(python_name<Float>("Matrix"),no_init);
    export_matrix_class<Float>(matrix_class);
    export_matrix_conversion<Float,Float>(matrix_class);
    export_matrix_arithmetic<Float,Float,Float>(matrix_class);
    export_matrix_arithmetic<Interval,Float,Interval>(matrix_class);
}

template<> void export_matrix<Interval>()
{
    class_< Matrix<Interval> > matrix_class(python_name<Interval>("Matrix"),no_init);
    export_matrix_class<Interval>(matrix_class);
    export_matrix_conversion<Interval,Float>(matrix_class);
    export_matrix_conversion<Interval,Interval>(matrix_class);
    export_matrix_arithmetic<Interval,Interval,Interval>(matrix_class);
    export_matrix_arithmetic<Interval,Interval,Float>(matrix_class);
    def("gs_inverse", (Matrix<Interval>(*)(const Matrix<Interval>&)) &gs_inverse);
    def("lu_inverse", (Matrix<Interval>(*)(const Matrix<Interval>&)) &lu_inverse);
    def("gs_solve", (Matrix<Interval>(*)(const Matrix<Interval>&,const Matrix<Interval>&)) &gs_solve);
    def("lu_solve", (Matrix<Interval>(*)(const Matrix<Interval>&,const Matrix<Interval>&)) &lu_solve);
    def("midpoint", (Matrix<Float>(*)(const Matrix<Interval>&)) &midpoint);
    implicitly_convertible< Matrix<Float>, Matrix<Interval> >();
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
