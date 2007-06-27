/***************************************************************************
 *            lpfsp.template.h
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
#ifndef ARIADNE_LPFSP_TEMPLATE_H
#define ARIADNE_LPFSP_TEMPLATE_H

namespace Ariadne {
  namespace LinearProgramming {
    
    template<class AP>
    bool lpfsp(const LinearAlgebra::Matrix<AP>& A, const LinearAlgebra::Vector<AP>& b,
    LinearAlgebra::Permutation& perm,
    LinearAlgebra::Vector<AP>& x, LinearAlgebra::Vector<AP>& y) {
      
      size_type m = A.number_of_rows();
      size_type n = A.number_of_columns();
      uint leave;
      
      if (x.size() != n)
        x = LinearAlgebra::Vector<AP>(n);
      
      // Check for feasibility of origin
      // if infeasibleorigin find most negative b; the related slack will leave the basis
      AP min_constr = 0;
      bool feasibleorigin = true;
      for (uint i = 0; i < b.size(); i++) // no break in loop to find most negative b
        if (b[i] < min_constr) {
          feasibleorigin = false;
          min_constr = b[i];
          leave = i;
        }
      
      if (feasibleorigin) {
        if (verbosity>2)
          std::clog << "Feasible Origin" << std::endl;
        for (uint i = 0; i < n; i++)
          x[i] = 0;
      }
      else { // solve the auxiliary problem
        if (verbosity>2)
          std::clog << "Will solve auxiliary problem" << std::endl;
        // add an auxiliary variable to the tableau; via a matrix operation
        LinearAlgebra::Matrix<AP> Aux = LinearAlgebra::Matrix<AP>(m, 1);
        LinearAlgebra::Permutation p = LinearAlgebra::Permutation(m+n+1);
        for (uint i = 0; i < m; i++)  Aux(i, 0) = -1;
        LinearAlgebra::Matrix<AP>R =  LinearAlgebra::concatenate_columns(Aux, A);
        uint auxcolumn = 0; // in the previous two lines we inserted an auxiliary variable in column 0
        LinearAlgebra::Vector<AP>c = LinearAlgebra::Vector<AP>(n+1);
        c(0) = -1;
        
        LinearAlgebra::Matrix<AP> t = to_tableau<AP>(R, b, c);
        uint rinc = t.row_increment();
        uint cinc = t.column_increment();
        LinearAlgebra::Matrix<AP> B = R;
        uint Brinc = B.row_increment();
        uint Bcinc = B.column_increment();
        AP* td = t.begin();
        AP* bptr = td + (m+n+1)*cinc;
        AP* cptr = td+m*rinc;
        AP* dptr = td + m*rinc + (m+n+1)*cinc;
        
        // auxiliary variable (x0) enters the basis, the slack of the most negative b leaves
        pivot_tableau(m, n+1, td, rinc, cinc, bptr, rinc, cptr, cinc, dptr, p, auxcolumn, leave);
        
        // perform simplex steps (note: the last parameter indicates it is an auxiliary problem)
        while (lpstp(m, n+1, td, rinc, cinc, bptr, rinc, cptr, cinc, dptr, p, B.begin(), Brinc, Bcinc, auxcolumn) == false);
        
        if (*dptr == 0) { // infeasible origin, but feasible solution
          // fill the x vector
          for (uint column = 0, i=0; column < n+1; column++)
            if (column != auxcolumn) {
              // find in which columns the decision variables can be found
              int j, index = p.getindex(column);
              // decision variable in base?
              if (index >= n+1) {
                for (j = 0; j < m; j++)
                  if (*(td + index*cinc + j*rinc) > 0) break;
                x[i] = *(bptr + j*rinc);
              }
              else
                x[i] = 0;
              i++;
            }
        }
        else // infeasible solution
          return false;
      }
      return true;
    }
    
    
    template<class AP>
    bool lpfsc(const LinearAlgebra::Matrix<AP>& A, const LinearAlgebra::Vector<AP>& b,
    LinearAlgebra::Vector<AP>& l, LinearAlgebra::Vector<AP>& u,
    LinearAlgebra::Permutation& perm,
    LinearAlgebra::Vector<AP>& x, LinearAlgebra::Vector<AP>& y) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class AP>
    bool lpfsd(const LinearAlgebra::Matrix<AP>& A, const LinearAlgebra::Vector<AP>& c,
    LinearAlgebra::Permutation& perm,
    LinearAlgebra::Vector<AP>& x, LinearAlgebra::Vector<AP>& y) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class R, class AP>
    tribool lprfsp(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b,
    LinearAlgebra::Permutation& perm,
    LinearAlgebra::Vector<AP>& x, LinearAlgebra::Vector<AP>& y) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class R, class AP>
    tribool lprfsc(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b,
    LinearAlgebra::Vector<R>& l, LinearAlgebra::Vector<R>& u,
    LinearAlgebra::Permutation& perm,
    LinearAlgebra::Vector<AP>& x, LinearAlgebra::Vector<AP>& y) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class R, class AP>
    tribool lprfsd(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& c,
    LinearAlgebra::Permutation& perm,
    LinearAlgebra::Vector<AP>& x, LinearAlgebra::Vector<AP>& y) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
  }//namespace LinearProgramming
}//namespace Ariadne

#endif /* ARIADNE_LPSFP_TEMPLATE_H */
