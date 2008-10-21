/***************************************************************************
 *            matrix.h
 *
 *  Copyright 2005-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file matrix.h
 *  \brief Matrices.
 */

#ifndef ARIADNE_MATRIX_H
#define ARIADNE_MATRIX_H 

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "numeric.h"

namespace Ariadne {

/// A matrix over a field.
template<class X>
class Matrix
  : public boost::numeric::ublas::matrix<X>
{
 public:
  Matrix()
    : boost::numeric::ublas::matrix<X>() { }
  template<class T> Matrix(const T& t)
    : boost::numeric::ublas::matrix<X>(t) { }
  Matrix(uint r, uint c)
    : boost::numeric::ublas::matrix<X>(r,c) { for(uint i=0; i!=r; ++i) { for(uint j=0; j!=c; ++j) { (*this)(i,j)=0; } } }
  template<class XX> Matrix(uint r, uint c, const XX* ptr)
    : boost::numeric::ublas::matrix<X>(r,c) { for(uint i=0; i!=r; ++i) { for(uint j=0; j!=c; ++j) { (*this)(i,j)=ptr[i*c+j]; } } }
  uint row_size() const { return this->size1(); }
  uint column_size() const { return this->size2(); }
  template<class XX> bool operator==(const Matrix<XX>& mx) const;
  const X* operator[](uint r) const { return &this->operator()(r,0); }
  X* operator[](uint r) { return &this->operator()(r,0); }
  const X& get(uint i, uint j) const { return (*this)[i][j]; }
  template<class T> void set(uint i, uint j, const T& x) { (*this)[i][j] = x; }
};

template<class X> Matrix<X> inverse(const Matrix<X>& A);
template<class X> X norm(const Matrix<X>& A);

template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& v) {
  return boost::numeric::ublas::prod(A,v);
}
 
inline Vector<Interval> operator*(const Matrix<Float>& A, const Vector<Interval>& v) {
  return boost::numeric::ublas::prod(A,v);
}

template<class X> inline Matrix<X> operator*(const Matrix<X>& A, const Matrix<X>& B) {
  return boost::numeric::ublas::prod(A,B);
}
 
template<class X> template<class XX> bool Matrix<X>::operator==(const Matrix<XX>& A2) const 
{
  const Matrix<X>& A1=*this;
  if(A1.row_size()!=A2.row_size() || A1.column_size() != A2.column_size()) {
    return false;
  }
  for(uint i=0; i!=A1.row_size(); ++i) {
   for(uint j=0; j!=A1.column_size(); ++j) {
     if(A1[i][j]!=A2[i][j]) {
       return false;
     }
   }
  }
  return true;
}



template<class X> Matrix<X> inverse(const Matrix<X>& A) 
{
  assert(A.row_size()==A.column_size());
  if(A.row_size()==1) {
    Matrix<X> I(1,1);
    I[0][0]=1/A[0][0];
    return I;
  }
  if(A.row_size()==2) {
    Matrix<X> I(2,2);
    X rdet = 1/(A[0][0]*A[1][1] - A[0][1] * A[1][0]);
    I[0][0]=A[1][1]*rdet;
    I[0][1]=-A[0][1]*rdet;
    I[1][0]=-A[1][0]*rdet;
    I[1][1]=A[0][0]*rdet;
    return I;
  }
  assert(false);
}

template<class X> X norm(const Matrix<X>& A) 
{
  X result=0;
  for(uint i=0; i!=A.row_size(); ++i) {
    X row_sum=0;
    for(uint j=0; j!=A.column_size(); ++j) {
      row_sum+=abs(A[i][j]);
    }
    if(row_sum>result) { result=row_sum; }
  }
  return result;
}

template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A) {
  if(A.row_size()==0 || A.column_size()==0) { os << '['; }
  for(uint i=0; i!=A.row_size(); ++i) { 
    for(uint j=0; j!=A.column_size(); ++j) { 
      os << (j==0 ? (i==0 ? '[' : ';') : ',') << A(i,j); } }
  return os << ']';
}


Matrix<Float> midpoint(const Matrix<Interval>& A);

} // namespace Ariadne

#endif
