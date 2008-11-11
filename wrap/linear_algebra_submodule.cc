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
 
#include "numeric.h"
#include "vector.h"
#include "matrix.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;

template<class X0, class X1, class X2>
Vector<X0> __mul__(const Matrix<X1>& A1, const Vector<X2>& v2) {
    return prod(A1,v2); 
}

template<class X0, class X1, class X2>
Matrix<X0> __mul__(const Matrix<X1>& A1, const Matrix<X2>& A2) {
    return prod(A1,A2); 
}



template<class X>
void export_vector()
{
    typedef Vector<X> V;

    class_< Vector<X> > vector_class("Vector",init<int>());
    //vector_class.def(init<const boost::python::object&>());
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
}


template<class X>
void export_matrix()
{
    class_< Matrix<X> > matrix_class("Matrix",init<size_t,size_t>());
    //matrix_class.def(init<const boost::python::object&>());
    matrix_class.def("rows", &Matrix<X>::row_size);
    matrix_class.def("columns", &Matrix<X>::column_size);
    matrix_class.def("row_size", &Matrix<X>::row_size);
    matrix_class.def("column_size", &Matrix<X>::column_size);
    matrix_class.def("__setitem__", (void(Matrix<X>::*)(size_t,size_t,const double&)) &Matrix<X>::set);
    matrix_class.def("__setitem__", (void(Matrix<X>::*)(size_t,size_t,const X&)) &Matrix<X>::set);
    matrix_class.def("__getitem__", &Matrix<X>::get, return_value_policy<copy_const_reference>());
    matrix_class.def(-self);        // __neg__
    matrix_class.def(self + self);  // __add__
    matrix_class.def(self - self);  // __sub__
    matrix_class.def(X() * self);  // __mul__
    //matrix_class.def(self * X());  // __mul__
    matrix_class.def(self / X());  // __div__
    matrix_class.def(double() * self);  // __mul__
    //matrix_class.def(self * double());  // __mul__
    matrix_class.def(self / double());  // __div__
    matrix_class.def("__mul__",(Vector<X>(*)(const Matrix<X>&,const Vector<X>&))&__mul__);
    matrix_class.def("__mul__",(Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&))&__mul__);
    //matrix_class.def("inverse", &inverse<X>);
    //matrix_class.def("determinant", &Matrix::determinant);
    //matrix_class.def("transpose", &Matrix::transpose);
    //matrix_class.def("solve", &Matrix::solve);
    matrix_class.def(boost::python::self_ns::str(self));    // __str__
}

template void export_vector<Float>();
template void export_vector<Interval>();
template void export_matrix<Float>();
template void export_matrix<Interval>();

void linear_algebra_submodule() {
    //export_vector<Float>();
    export_vector<Interval>();
    //export_matrix<Float>();
    export_matrix<Interval>();
}
