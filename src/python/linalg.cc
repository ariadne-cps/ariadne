/***************************************************************************
 *            python/linalg.cc
 *
 *  21 October 2005
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
 */// g++ -fpic -I../include/ -I/usr/include/python2.4 -c interval.cc
// g++ -shared interval.o -o interval.so -lboost_python -lmpfr -lgmp  -lsuperlu -lblas

#include <iostream>

#include <boost/python.hpp>

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

#include <sla/sla.h>

typedef SLA::Vector<Interval> Vector;
typedef SLA::Matrix<Interval> Matrix;

typedef Vector vector_unary_func(const Vector& v);
typedef Vector vector_binary_func(const Vector& v, const Vector& v);

Vector operator-(const Vector& v) {
  return SLA::operator-(v); }
Vector operator+(const Vector& v1, const Vector& v2) {
  return SLA::operator+(v1,v2); }
Vector operator-(const Vector& v1, const Vector& v2) {
  return SLA::operator-(v1,v2); }
Vector operator*(const Interval& i, const Vector& v) {
  return SLA::operator*(v,i); }
Vector operator*(const Vector& v, const Interval& i) {
  return SLA::operator*(v,i); }
Vector operator/(const Vector& v, const Interval& i) {
  return SLA::operator/(v,i); }

Matrix operator-(const Matrix& a) {
  return SLA::operator-(a); }
Matrix operator+(const Matrix& a1, const Matrix& a2) {
  return SLA::operator+(a1,a2); }
Matrix operator-(const Matrix& a1, const Matrix& a2) {
  return SLA::operator-(a1,a2); }
Matrix operator*(const Interval& i, const Matrix& a) {
  return SLA::operator*(a,i); }
Matrix operator*(const Matrix& a, const Interval& i) {
  return SLA::operator*(a,i); }
Vector operator*(const Matrix& a, const Vector& v) {
  return SLA::operator*(a,v); }
Matrix operator*(const Matrix& a, const Matrix& v) {
  return SLA::operator*(a,v); }
Matrix operator/(const Matrix& a, const Interval& i) {
  return SLA::operator/(a,i); }


BOOST_PYTHON_MODULE(linalg)
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    typedef Interval (*IFUN)(const Interval&);

    class_< Vector >("Vector")
      .def(init<int>())
      .def(-self)        // __neg__
      .def(self + self)  // __add__
      .def(self - self)  // __sub__
      .def(Interval() * self)  // __mul__
      .def(self * Interval())  // __mul__
      .def(self / Interval())  // __div__
      .def(boost::python::self_ns::str(self))    // __str__
      .def("__setitem__", &Vector::set)
      .def("__getitem__", &Vector::get, return_value_policy<copy_const_reference>())
      ;

    class_< Matrix >("Matrix")
      .def(init<int,int>())
      .def(-self)        // __neg__
      .def(self + self)  // __add__
      .def(self - self)  // __sub__
      .def(Interval() * self)  // __mul__
      .def(self * Interval())  // __mul__
      .def(self * Vector())    // __mul__
      .def(self * Matrix())    // __mul__
      .def(self / Interval())  // __div__
      .def(boost::python::self_ns::str(self))    // __str__
      .def("set", &Matrix::set)
      .def("get", &Matrix::get, return_value_policy<copy_const_reference>())
      .def("inverse", &Matrix::inverse)
      .def("determinant", &Matrix::determinant)
      .def("transpose", &Matrix::transpose)
      .def("solve", &Matrix::solve)
      .def("__setitem__", &Matrix::set)
      .def("__getitem__", &Matrix::get, return_value_policy<copy_const_reference>())
      ;


}
