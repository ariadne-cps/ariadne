/***************************************************************************
 *            lpslv.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
 ****************************************************************************/
/*
 * Based on the linear programming algorithms in PPL-0.8
 *   Copyright (C) 2001-2006 Roberto Bagnara <bagnara@cs.unipr.it>
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
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file lpslv.h
 *  \brief Linear programming solver.
 */

#ifndef _ARIADNE_LPSLV_H
#define _ARIADNE_LPSLV_H

#include <iosfwd>
#include <cassert>
#include <map>


#include "../linear_algebra/matrix.h"

namespace Ariadne { namespace LinearAlgebra {

int verbosity=0;
  
/*! \ingroup LinearAlgebra
 *  \brief Solver for linear programming problems.
 *
 *  \param m Number of free constraints
 *  \param n Number of free variables
 *  \param (A,rincA,cincA) An m-by-n matrix
 *  \param (B,incB) An m-element vector
 *  \param (C,incC) An n-element vector
 *  \param d A scalar
 *  \param piv A consecutive (m+n)-element vector of integers.
 *
 * Solve the linear programming problem 
 * \f$ text{minimize}\ c^Tx \text{ subject to } Ax+y=b,\ x,y\geq0\f$ where the
 * current point is given by \f$x=0,\ y=b\f$ and the current value is \f$-d\f$. 
 */
template<class R>
void lpslv(int m, int n, R* A, int rincA, int cincA, R* B, int incB, R* C, int incC, R& d, int* piv)
{
  if(verbosity>1) {
    std::cerr << "lpslv(" << m << "," << n << ", " << A-A << "," << rincA << "," << cincA << ", " 
              << B-A << "," << incB << ", " << C-A << "," << incC << ", "
              << &d-A << ", " << piv << ")" << std::endl;
    std::cerr << "T=" << Matrix<R>(m+1,n+1,A,rincA,cincA) << std::endl;
    std::cerr << "A=" << Matrix<R>(m,n,A,rincA,cincA) << "; b=" << Vector<R>(m,B,incB)
              << "; c=" << Vector<R>(n,C,incC) << "; d=" << d << std::endl;
  }
  
  R one=static_cast<R>(1);
  size_type recursions=16;
  int i,j;
      
  //Select variable to enter basis
  for(j=0; j!=n; ++j) {
    if(C[j*incC] < 0) {
      break;
    }
  }

  while(j!=n) {
    assert(--recursions);
    
    // compute variable to exit basis
    i=m;
    R min_change=0;
    for(int k=0; k!=m; ++k) {
      if(A[k*rincA+j*cincA]>0) {
        R change=B[k*incB]/A[k*rincA+j*cincA];
        if(change<min_change || min_change==0) {
          min_change=change;
          i=k;
        }
      }
    }

    if(verbosity>1) {
      std::cerr << "Pivoting on (exit=" << i << ", enter=" << j << ")" << std::endl;
    }
    
    std::swap(piv[j],piv[n+i]);
      
    // Modify the tableau
    R pivot=A[i*rincA+j*cincA];
    
    // Subtract A(p,j)/A(i,j) times row i from row p, p!=i,
    // except in the jth column, which is divided by A(i,j)
    for(int p=0; p!=m; ++p) {
      if(p!=i) {
        R scale=A[p*rincA+j*cincA]/pivot;
        for(int q=0; q!=n; ++q) {
          if(q!=j) {
            A[p*rincA+q*cincA] -= A[i*rincA+q*cincA]*scale;
          }
        }
        B[p*incB] -= B[i*incB]*scale;
        A[p*rincA+j*cincA] = -scale;
      }
    }
    // Subtract c(j)/A(i,j) times row i from c,
    // except in the jth column, which is divided by A(i,j)
    R scale = C[j*incC]/pivot;
    for(int q=0; q!=n; ++q) {
      if(q!=j) {
        C[q*incC] -= A[i*rincA+q*incC]*scale;
      }
    }
    d -= B[i*incB]*scale;
    C[j*incC] = -scale;

    // Scale the ith row by 1/pivot, except the jth column, which is set to 1/pivot
    scale=one/pivot;
    for(int q=0; q!=n; ++q) {
      A[i*rincA+q*cincA] *= scale;
    }
    B[i*incB] *= scale;
    A[i*rincA+j*cincA] = scale;
    
    if(verbosity>1) {
      std::cerr << "A=" << Matrix<R>(m,n,A,rincA,cincA) << "; b=" << Vector<R>(m,B,incB) 
                << "; c=" << Vector<R>(n,C,incC) << "; d=" << d << std::endl;
    }
    
    // Select variable to enter basis
    for(j=0; j!=n; ++j) {
      if(C[j*incC] < 0) {
        break;
      }
    }
  }
  
  return;

}

}}


#endif /* _ARIADNE_LPSLV_H */
