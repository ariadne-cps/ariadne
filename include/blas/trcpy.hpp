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

#ifndef __BLAS_TRCPY_HPP__
#define __BLAS_TRCPY_HPP__

#include "blas.hpp"
 
template<typename scalarA, typename scalarB>
void
BLAS::trcpy(ORDER order, UPLO uplo, DIAG diag, int m, int n, 
            const scalarA *A, int ldA, scalarB *B, int ldB)
{
  if ( (order==RowMajor && uplo==Lower) || (order==ColMajor && uplo==Upper)) {
    for(int i=0; i!=m; ++i) {
      for(int j=0; j!=i; ++j) {
        B[i*ldB+j]=A[i*ldA+j];
      }
    }
  } else {  
    for(int i=0; i!=m; ++i) {
      for(int j=i+1; j!=n; ++j) {
        B[i*ldB+j]=A[i*ldA+j];
      }
    }
  }
  
  if(diag==Unit) {
    for(int i=0; i!=min(m,n); ++i) {
      B[i*ldB+i]=1;
    }
  } else {
    for(int i=0; i!=min(m,n); ++i) {
      B[i*ldB+i]=A[i*ldA+i];
    }
  }
  
  return;
}

#endif // __BLAS_TRCPY_HPP__
