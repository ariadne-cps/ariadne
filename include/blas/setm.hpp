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

#ifndef __BLAS_SETM_HPP__
#define __BLAS_SETM_HPP__

#include "blas.hpp"

template<typename scalar>
void
BLAS::setm (const enum ORDER Order, const int M, const int N, 
            const scalar alpha, scalar *A, const int ldA)
{

  int i,j;
  if(Order==RowMajor) {
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        A[ldA*i+j] = alpha;
      }
    }
  } else {
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        A[i+ldA*j] = alpha;
      }
    }
  }

}

  

#endif // __BLAS_SETM_HPP__
