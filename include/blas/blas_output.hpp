/*
 * Copyright (C) 2006 Pieter Collins <Pieter.Collins@cwi.nl>
 *
 * Based on the BLAS implementation in Gnu Scientific Library 1.8
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 */

/*  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*! \file blas_output.hh
 *  \brief Output routines for BLAS vectors and matrices.
 */

#ifndef __BLAS_OUTPUT_HPP__
#define __BLAS_OUTPUT_HPP__

#include <iostream>
#include "blas.hpp"

namespace BLAS {

/*
 * ===========================================================================
 * Output routines
 * ===========================================================================
 */

template<typename scalar> class Vector {
 public:
  Vector(const int N, const scalar *X, const int incX)
    : n(N), x(X), incx(incX) { }
  int size() const { return n; }
  const scalar& operator[] (const int& i) const { return x[i*incx]; }
 private:
  const int n; const scalar *x; const int incx;
};

template<typename scalar> class Matrix {
 public:
  Matrix(const int &M, const int &N, const scalar *A, const int &ldA)
    : m(M), n(N), a(A), lda(ldA) { }
  int number_of_rows() const { return m; }
  int number_of_columns() const { return n; }
  const scalar* operator[] (const int& i) const { return &a[i*lda]; }
 private:
  const int m; const int n; const scalar *a; const int lda;
};

template<typename scalar> inline 
Vector<scalar> 
vector(const int N, const scalar *X, const int incX) { 
  return Vector<scalar>(N,X,incX); }

template<typename scalar> inline 
Vector<scalar> 
vector(const int N, const scalar *X) {
  return Vector<scalar>(N,X,1); }

template<typename scalar> inline 
Matrix<scalar> 
matrix(const int M, const int N, const scalar *A, const int ldA) { 
  return Matrix<scalar>(M,N,A,ldA); }

template<typename scalar> inline 
Matrix<scalar> 
matrix(const int M, const int N, const scalar *A) { 
  return Matrix<scalar>(M,N,A,N); }

template<typename scalar> inline 
std::ostream& 
operator<< (std::ostream& os, const Vector<scalar>& v) 
{
  for(int i=0; i!=v.size(); ++i) {
    os << ( i==0 ? "[" : "," ) << v[i];
  }
  os << "]";
  return os;
}

template<typename scalar> inline 
std::ostream& 
operator<< (std::ostream& os, const Matrix<scalar>& A) 
{
  for(int i=0; i!=A.number_of_rows(); ++i) {
    for(int j=0; j!=A.number_of_columns(); ++j) {
      os << ( j==0 ? "[" : "," ) << A[i][j];
    }
    os << "]\n";
  }
  return os;
}

template<typename real> inline 
std::ostream& 
operator<< (std::ostream& os, const complex<real>& z) 
{
  if(z.real()!=0) { os << z.real(); if(z.imag()>0) { os << "+"; } }
  if(z.imag()!=0) { os << z.imag() << "i"; }
  return os;
}

}

#endif /* __BLAS_OUTPUT_HPP__ */

