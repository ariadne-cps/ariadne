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
#include <map>


#include "../linear_algebra/matrix.h"

namespace Ariadne { namespace LinearAlgebra {

/*! \ingroup LinearAlgebra
 *  \brief Solver for linear programming problems.
 *
 *  \param m Number of free constraints
 *  \param n Number of free variables
 */
template<typename R>
void lpslv(int m, int n, R* A, int rincA, int cincA, R* b, int incB, R* c, int incC, R& d, int* piv)
{
  R one=static_cast<R>(1);
  size_type recursions=100;
  int i,j;
      
  //Select variable to enter basis
  for(j=0; j!=n; ++j) {
    if(c[j*incC] < 0) {
      break;
    }
  }

  while(j!=n) {
    assert(--recursions);
    
    // compute variable to exit basis
    i=m;
    R min_change=0;
    for(size_type k=0; k!=m; ++k) {
      if(A[k*rincA+j*cincA]>0) {
        R change=b[k*incB]/A[k*rincA+j*cincA];
        if(change<min_change || min_change==0) {
          min_change=change;
          i=k;
        }
      }
    }

    // std::cerr << "Pivoting on (" << i << "," << j << ")" << std::endl;
        
    std::swap(piv[j],piv[n+i]);
      
    // Modify the tableau
    R pivot=A[i*rincA+j*cincA];
    
    // Subtract A(p,j)/A(i,j) times row i from row p,
    // except in the jth column, which is divided by A(i,j)
    for(int p=0; p!=m; ++p) {
      if(p!=i) {
        R scale=A[p*rincA+j*cincA]/pivot;
        for(int q=0; q!=n; ++q) {
          A[p*rincA+q*cincA] -= A[i*rincA+q*cincA]*scale;
        }
        b[p*incB] -= b[i*incB]*scale;
        A[p*rincA+j*cincA] = scale;
      }
    }
    // Subtract c(j)/A(i,j) times row i from c,
    // except in the jth column, which is divided by A(i,j)
    R scale = c[j*incC]/pivot;
    for(size_type q=0; q!=n; ++q) {
      c[q*incC] -= c[j*incC]*scale;
    }
    d -= c[j*incC]*scale;
    c[j*incC] = scale;

    // Set the ith row
    scale=one/pivot;
    for(size_type q=0; q!=n; ++q) {
      A[i*rincA+q*cincA] *= scale;
    }
    b[i*incB] *= scale;
    A[i*rincA+j*cincA] = scale;
    
    // Select variable to enter basis
    for(j=0; j!=n; ++j) {
      if(c[j*incC] < 0) {
        break;
      }
    }
  }
  
  return;

}

}}


#endif /* _ARIADNE_LPSLV_H */
