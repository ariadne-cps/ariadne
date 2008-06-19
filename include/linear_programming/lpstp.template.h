/***************************************************************************
 *            lpstp.template.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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

#include "lputil.h"

namespace Ariadne {





  template<class AP>
  tribool 
  lpstp(const Matrix<AP>& A, 
                           const Vector<AP>& b, 
                           const Vector<AP>& c,
                           Permutation& p, 
                           Matrix<AP>& B,
                           // Use these variables too if it's more efficient
                           Vector<AP>& x, 
                           Vector<AP>& y, 
                           Vector<AP>& z) 
  {
      
    if (verbosity > 2)
      std::clog << "A:" << A << ", b:" << b << ", c:" << c << ", p:" << p << ", B:" << B << std::endl;
      
    Matrix<AP> t = to_tableau<AP>(A, b, c);
      
    uint m = A.number_of_rows();
    uint n = A.number_of_columns();
    uint rinc = t.row_increment();
    uint cinc = t.column_increment();
    
    if (verbosity > 3)
      std::clog << "t:" << t << std::endl << "m:" << m << ", n:" << n << ", ri:" << rinc << ", ci:" << cinc << std::endl;
    
    uint Brinc = B.row_increment();
    uint Bcinc = B.column_increment();
    
    AP* td = t.begin();
    uint* pptr = p.begin();
    return lpstp(m, n, td, rinc, cinc, td + (n+m)*cinc, rinc, td + m*rinc, cinc, td + m*rinc + (n+m)*cinc, pptr, B.begin(), Brinc, Bcinc);
  }
    
    
  template<class AP>
  tribool 
  lpstpc(const Matrix<AP>& A, 
                            const Vector<AP>& b, 
                            const Vector<AP>& c,
                            const Vector<AP>& l, 
                            const Vector<AP>& u,
                            Permutation& p, 
                            Matrix<AP>& B, 
                            Vector<AP>& x) 
  {
    throw NotImplemented(__PRETTY_FUNCTION__);
  }
    
    

  template<class AP>
  tribool 
  lpstp(uint m, uint n,
                           AP* Aptr, uint Arinc, uint Acinc,
                           AP* bptr, uint binc,
                           AP* cptr, uint cinc,
                           AP* dptr,
                           uint* pptr,
                           AP* Bptr, uint Brinc, uint Bcinc
                           // use only if it's more efficient
                           // , AP* xptr, uint xinc, AP* yptr, uint yinc, AP* zptr, uint zinc
                           ) 
  {
    return lpstp(m, n, Aptr, Arinc, Acinc, bptr, binc, cptr, cinc, dptr, pptr, Bptr, Brinc, Bcinc, -1);
  }
    

  template<class AP>
  tribool 
  lpstp(uint m, uint n,
                           AP* Aptr, uint Arinc, uint Acinc,
                           AP* bptr, uint binc,
                           AP* cptr, uint cinc,
                           AP* dptr,
                           uint* pptr,
                           AP* Bptr, uint Brinc, uint Bcinc,
                           // use only if it's more efficient
                           // AP* xptr, uint xinc, AP* yptr, uint yinc, AP* zptr, uint zinc,
                           int aux // index of additional variable to solve the auxiliary problem
                          )
  {
    if (verbosity > 2) {
      std::clog << "lpstp(m=" << m << ", n=" << n << ", Arinc=" << Arinc << ", Acinc=" << Acinc
                << ", binc=" << binc << ", cinc=" << cinc << ", perm=" << Permutation(pptr,pptr+n) << ", aux=" << aux << ")" << std::endl;
      std::clog << "tableau=" << Matrix<AP>(m+1, m+n+1, Aptr, Arinc, Acinc) << std::endl;
    }
      
    int leave, enter = -1;
    
    // select variable to enter basis, strategy for auxiliary problem differs
    if (aux >= 0) {
      // auxiliary problem -> pick the one having the most negative coefficient
      AP minval = 0;
      for (uint i = 0; i < n; i++) {
        AP val = cptr[i*cinc];
        if (val < minval) {
          minval = val;
          enter = i;
        }
      }
    } else {
      // pick the first one having a negative coefficient
      for (uint i = 0; i < n; i++)
        if (cptr[i*cinc] < 0) {
          enter = i;
          break;
        }
    }
      
    if (enter < 0) {
      if (verbosity > 2) std::clog << "lpstp: no variable to enter basis -- current solution is optimal" << std::endl;
      return true;  // current solution is optimal
    }
    
    // find the row of the auxiliary variable (if present) to make sure we can get it out of the basis
    int auxrow = -1;
    if (aux >= 0) {
      leave = getindex(pptr,n,aux);
      if (leave > n) {  // auxiliary variable is indeed part of the basis
        for (uint row = 0; row < m; row++)
          if (Aptr[leave*Acinc + row*Arinc] > 0) {
            auxrow = row;
            break;
          }
      }
    }
    
    // compute variable to leave basis; i.e. the one having the most stringent upper bound on increase of the entering variable
    leave = m;
    AP min_ratio = -1; // initialize: note that ratio always >= 0
    for (uint k = 0; k < m; k++) {
      int indx = k*Arinc + enter*Acinc;
      if (Aptr[indx] > 0) {
        AP ratio = bptr[k*binc] / Aptr[indx];
        if (ratio < min_ratio || min_ratio == -1 || (ratio == min_ratio && k == auxrow)) {
          min_ratio = ratio;
          leave = k;
        }
      }
    }
    if (min_ratio == -1) {
      if (verbosity > 1) std::clog << "lpstp: no positive entry in pivot column -- problem is unbounded" << std::endl;
      return indeterminate; // problem is unbounded
    }
    pivot_tableau(m, n, Aptr, Arinc, Acinc, bptr, binc, cptr, cinc, dptr, pptr, enter, leave);
    return false;
  }
    
        

  template<class AP>
  tribool 
  lpstpc(uint m, uint n,
                            AP* Aptr, uint Arinc, uint Acinc,
                            AP* bptr, uint binc,
                            AP* cptr, uint cinc,
                            AP* lptr, uint linc,
                            AP* uptr, uint uinc,
                            AP* dptr,
                            uint* pptr,
                            AP* Bptr, uint Brinc, uint Bcinc
                            // , AP* xptr, uint xinc, AP* yptr, uint yinc, AP* zptr, uint zinc
                            ) 
  {
    throw NotImplemented(__PRETTY_FUNCTION__);
  }
    



} // namespace Ariadne
                                                                                                                                  
