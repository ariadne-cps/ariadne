/***************************************************************************
 *      lpslv.template.h
 *
 * Copyright 2006 Alberto Casagrande, Pieter Collins
 * casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <cassert>
#include "linear_algebra/exceptions.h"

#include "lputil.h"
#include "lpfsp.h"
#include "lpstp.h"

namespace Ariadne {
  
    
    
    template<class R, class AP>
    AP lpslv(const Matrix<R>& A, 
             const Vector<R>& b, 
             const Vector<R>& c,
             Permutation& p, 
             Vector<AP>& x, 
             Vector<AP>& y) 
    {
      size_type m = A.number_of_rows();
      size_type n = A.number_of_columns();
      Matrix<R>t;
      Matrix<R>B;
      uint rinc, cinc;
      uint Brinc, Bcinc;
      
      R* td;
      R* bptr;
      R* cptr;
      R* dptr;
      uint* pptr;
      tribool ans;
      
      // check feasibility
      if (!lpfsp(A, b, p, x, y))
        ARIADNE_THROW(InfeasibleSolution, __PRETTY_FUNCTION__, "lp problem is infeasible");
      
      // initialize variables
      if (x == Vector<R>::zero(n)) {        // setup for feasible origin
        t = to_tableau<R>(A, b, c);
        rinc = t.row_increment();
        cinc = t.column_increment();
        B = A;                                        // currently B is not used
        Brinc = B.row_increment();
        Bcinc = B.column_increment();
        td = t.begin();
        bptr = td + (n+m)*cinc;
        cptr = td+m*rinc;
        dptr = td + m*rinc + (m+n)*cinc;
        pptr = p.begin();
      }
      else {        // setup for infeasible origin
        
        // TODO
        ARIADNE_THROW(InfeasibleOrigin, __PRETTY_FUNCTION__, "lp problem with infeasible origin not implemented yet");
        
      }
      
      
      // simplex steps
      do {
        ans = lpstp(m, n, td, rinc, cinc, bptr, rinc, cptr, cinc, dptr, pptr, B.begin(), Brinc, Bcinc);
        // ToDo: check if d increments and if not then limit recursion
      } while (ans == false);
      
      if (x.size() != c.size())
        x = Vector<AP>(c.size());
      
      if (ans == true) { // found optimal solution
        
        // create x vector
        for (uint i = 0; i < x.size(); i++) {
          // find in which columns the decision variables can be found
          uint j;
          int index = p.getindex(i);
          // decision variable in base?
          if (index >= x.size()) {
            for (j = 0; j < b.size(); j++)
              if (*(td + index*cinc + j*rinc) > 0) break;
            x[i] = *(bptr + j*rinc);
          }
          else
            x[i] = 0;
        }
        
        // create y vector
        if (y.size() != b.size())
          y = Vector<AP>(b.size());
        for (int i = 0; i < y.size(); i++) {
          // find in which column the slack variables can be found
          int index = p.getindex(i + c.size());
          y[i] = t(m, index);
        }
        
      }
      else { // ans == indeterminate -> problem is unbounded
        // ToDo: create feasible x
        y = Vector<AP>();
      }
      
      if (verbosity > 2)
        std::clog << "d:" << *dptr << ", x:" << x << "\n" << std::endl;
      return *dptr;
    }
    
    

    template<class R, class AP>
    AP 
    lpslvc(const Matrix<R>& A, 
           const Vector<R>& b, 
           const Vector<R>& c,
           const Vector<R>& l, 
           const Vector<R>& u,
           Permutation& p, 
           Vector<AP>& x, 
           Vector<AP>& y) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class R>
    void lpslv(int m, int n, 
               R* A, int rincA, int cincA, 
               R* B, int incB, 
               R* C, int incC, 
               R& d, 
               int* piv) 
    {
      if (verbosity>1) {
        std::clog << "lpslv(" << m << "," << n << ", " << A-A << "," << rincA << "," << cincA << ", "
        << B-A << "," << incB << ", " << C-A << "," << incC << ", "
        << &d-A << ", " << piv << ")" << std::endl;
        std::clog << "T=" << Matrix<R>(m+1, n+1, A, rincA, cincA) << std::endl;
        std::clog << "A=" << Matrix<R>(m, n, A, rincA, cincA) << "; b=" << Vector<R>(m, B, incB)
        << "; c=" << Vector<R>(n, C, incC) << "; d=" << d << std::endl;
      }
      
      R one=static_cast<R>(1);
      size_type recursions=16;
      int i, j;
      
      //Select variable to enter basis
      for (j=0; j!=n; ++j) {
        if (C[j*incC] < 0) {
          break;
        }
      }
      
      while(j!=n) {
        --recursions;
        if(recursions==0) {
          std::cerr << "Error in lpslv(...) with A="<<Matrix<R>(m,n,A,rincA,cincA)<<", b="<<Vector<R>(m,B,incB)<<", c="<<Vector<R>(n,C,incC)<<std::endl;
          throw std::runtime_error("Unknown error in lpslv");
        }
        
        // compute variable to exit basis
        i=m;
        R min_change=-1;
        for (int k=0; k!=m; ++k) {
          if (A[k*rincA+j*cincA]>0) {
            R change=B[k*incB]/A[k*rincA+j*cincA];
            if (change<min_change || min_change== -1) {
              min_change=change;
              i=k;
            }
          }
        }
        
        if (verbosity>1) {
          std::clog << "Pivoting on (exit=" << i << ", enter=" << j << ")" << std::endl;
        }
        
        std::swap(piv[j], piv[n+i]);
        
        // Modify the tableau
        R pivot=A[i*rincA+j*cincA];
        
        // Subtract A(p,j)/A(i,j) times row i from row p, p!=i,
        // except in the jth column, which is divided by A(i,j)
        for (int p=0; p!=m; ++p) {
          if (p!=i) {
            R scale=A[p*rincA+j*cincA]/pivot;
            for (int q=0; q!=n; ++q) {
              if (q!=j) {
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
        for (int q=0; q!=n; ++q) {
          if (q!=j) {
            C[q*incC] -= A[i*rincA+q*incC]*scale;
          }
        }
        d -= B[i*incB]*scale;
        C[j*incC] = -scale;
        
        // Scale the ith row by 1/pivot, except the jth column, which is set to 1/pivot
        scale=one/pivot;
        for (int q=0; q!=n; ++q) {
          A[i*rincA+q*cincA] *= scale;
        }
        B[i*incB] *= scale;
        A[i*rincA+j*cincA] = scale;
        
        if (verbosity>1) {
          std::clog << "A=" << Matrix<R>(m, n, A, rincA, cincA) << "; b=" << Vector<R>(m, B, incB)
          << "; c=" << Vector<R>(n, C, incC) << "; d=" << d << std::endl;
        }
        
        // Select variable to enter basis
        for (j=0; j!=n; ++j) {
          if (C[j*incC] < 0) {
            break;
          }
        }
      }
      
      return;
    }
       

} // namespace Ariadne

