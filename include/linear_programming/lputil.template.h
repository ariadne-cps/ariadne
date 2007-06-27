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


template<class R>
LinearAlgebra::Matrix<R>
LinearProgramming::to_tableau(const LinearAlgebra::Matrix<R>& A, 
                              const LinearAlgebra::Vector<R>& b, 
                              const LinearAlgebra::Vector<R>& c) 
{
  size_type nc = b.size(); // number of constraints
  size_type nv = c.size(); // number of variables
  
  ARIADNE_CHECK_SIZE(b, A.number_of_rows(), __PRETTY_FUNCTION__)
    ARIADNE_CHECK_SIZE(c, A.number_of_columns(), __PRETTY_FUNCTION__)
    
    if (verbosity > 2)
      std::clog << "to_tableau in:" << std::endl << " A=" << A << std::endl << " b=" << b << std::endl << " c=" << c << std::endl;
  
  size_type n = nv+nc;
  // tableau  needs a row for each constraint + 1 and a column for each variable + slacks + 1
  LinearAlgebra::Matrix<R> tableau =  LinearAlgebra::Matrix<R>(nc+1, n+1);
  //  NB: uninitialized elements are assumed to be 0
  
  for (size_type i = 0; i != nc; i++) {
    
        // paste row of input matrix A
    for (size_type j = 0; j != nv; j++)
      tableau(i, j) = A(i, j);
    
    // slacks
    tableau(i, i + nv) = 1;
    
    // constraints
    tableau(i, n) = b(i);
  }
  
  // objective function
  for (size_type j = 0; j != nv; j++)
    tableau(nc, j) = -c(j);
  
  if (verbosity > 2)
    std::clog << "to_tableau out:" << tableau << std::endl;
  
  return tableau;
}



// Modify the tableau
template<class AP>
void 
LinearProgramming::pivot_tableau(uint m, uint n,
                                 AP* Aptr, uint Arinc, uint Acinc,
                                 AP* bptr, uint binc,
                                 AP* cptr, uint cinc,
                                 AP* dptr,
                                 uint* pptr,
                                 int enter, int leave) 
{
  
  int leave_inc = leave*Arinc;
  int enter_inc = enter*Acinc;
  
  AP pivot_scale = static_cast<AP>(1) / Aptr[leave_inc+enter_inc];
  
  if (verbosity > 3)
    std::clog << "pivot_tableau: leave=" << n+leave << ", enter=" << enter << ", pivot scale=" << pivot_scale << ", perm=" << LinearAlgebra::Permutation(pptr,pptr+n) << std::endl;
  
  std::swap(pptr[enter], pptr[n+leave]);
  
  // Subtract Aptr(p,enter) / Aptr(leave,enter) times row leave from row p, p!=leave,
  // except in the leave column, which is divided by Aptr(leave,enter)
  for (uint p=0; p < m; p++) {
    if (p != leave ) {
      uint p_inc = p*Arinc;
      AP scale = Aptr[p_inc + enter_inc] * pivot_scale;
      if (verbosity > 4)
        std::clog << "scale=" << scale << std::endl;
      for (uint q=0; q < n; q++) {
        if (q != enter ) {
          uint q_inc = q*Acinc;
          Aptr[p_inc + q_inc] -= Aptr[leave_inc + q_inc]*scale;
        }
      }
      bptr[p*binc] -= bptr[leave*binc]*scale;
      Aptr[p_inc + enter_inc] = -scale;
    }
  }
  
  if (verbosity > 4)
    std::clog << "after subtracting row " << leave << ": \ntableau=" << LinearAlgebra::Matrix<AP>(m+1, m+n+1, Aptr, Arinc, Acinc) << std::endl;
  
  // Subtract c(enter)/Aptr(leave,enter) times row leave from c,
  // except in the leave column, which is divided by Aptr(leave,enter)
  AP scale = cptr[enter*cinc] * pivot_scale;
  if (verbosity > 4)
    std::clog << "scale=" << scale << std::endl;
  for (uint q=0; q < n; q++) {
    if (q != enter) {
      uint q_inc = q*cinc;
      cptr[q_inc] -= Aptr[leave_inc+q_inc]*scale;
    }
  }
  *dptr -= bptr[leave*binc]*scale;
  cptr[enter*cinc] = -scale;
  
  if (verbosity > 4)
    std::clog << "after subtracting row " << leave << " from c: \ntableau=" << LinearAlgebra::Matrix<AP>(m+1, m+n+1, Aptr, Arinc, Acinc) << std::endl;
  
  // Scale the enter row, except the leave column
  scale = pivot_scale;
  for (uint q=0; q < n; q++) {
    Aptr[leave_inc+q*Acinc] *= scale;
  }
  bptr[leave*binc] *= scale;
  Aptr[leave_inc+enter_inc] = scale;
  
  if (verbosity > 2)
    std::clog << "leaving pivot_tableau\ntableau=" << LinearAlgebra::Matrix<AP>(m+1, m+n+1, Aptr, Arinc, Acinc) << "\nperm=" << LinearAlgebra::Permutation(pptr,pptr+n) << std::endl;
}


} // namespace Ariadne

