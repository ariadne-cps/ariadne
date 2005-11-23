/***************************************************************************
 *            linear_algebra.h
 *
 *  Mon May  3 12:31:15 2004
 *  Copyright  2004  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
/*! \file linear_algebra.h
 *  \brief Basic linear algebra.
 */

#ifndef _ARIADNE_LINEAR_ALGEBRA_H
#define _ARIADNE_LINEAR_ALGEBRA_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "numerical_type.h"

namespace Ariadne {
  namespace LinearAlgebra {

    using boost::numeric::ublas::vector;
    using boost::numeric::ublas::matrix;
  
    typedef size_t size_type;
        
    template <typename Real>
    inline vector<Real> null_vector(size_type dim) {
      vector<Real> v(dim);
      for (size_type j=0; j< dim; j++) {
        v(j)=0.0;
      }
      return v;
    }
    
    template <typename Real>
    inline matrix<Real> zero_matrix(size_type dim) {
      matrix<Real> A(dim,dim);
      for (size_type i=0; i<dim; ++i) {
        for (size_type j=0; j<dim; ++j) {
          A(i,j)=0.0;
        }
      }
      return A;
    }
    
    template <typename Real>
    inline matrix<Real> identity_matrix(size_type dim) {
      matrix<Real> A(dim,dim);
      for (size_type i=0; i<dim; i++) {
        for (size_type j=0; j<dim; j++) {
          A(i,j)=0.0;
        }
        A(i,i)=1.0;
      }
      return A;
    }
      
    template <typename Real>
    inline matrix<Real> exp_Ah(const matrix<Real> &A, 
                               const Real h, const unsigned int n) 
    {
      matrix<Real> tmp,e_Ah;
      
      e_Ah=identity_matrix<Real>(A.size1());
      tmp=e_Ah;
      
      /* tmp = \frac{h^{0}}{0!}*A^{0} = I
       * and e_Ah = \Sum_{j=0}^{0}\frac{h^j}{j!}*A^{j} = I */
      for (size_type i=1; i<n; ++i) {
        /* tmp = \frac{h^{i-1}}{(i-1)!}*A^{i-1}
         * and e_Ah = \Sum_{j=0}^{i-1}\frac{h^j}{j!}*A^{j} */
        tmp *= (h/i);
        tmp = prod(tmp,A);
        /* tmp =  (h^i/i!)*A^i */
        e_Ah += tmp;
        /*  e_Ah = \Sum_{j=0}^{i}\frac{h^j}{j!}*A^{j} */
      }
      return e_Ah;
    }
    
    template <typename Real> 
    inline vector<Real> exp_b_approx(const matrix<Real> &A, 
                                     const vector<Real> &b, 
                                     const Real h, const unsigned int n) 
    {
      matrix<Real> tmp,e_b;
      
      tmp=h*identity_matrix<Real>(A.size1());
      e_b=tmp;
      /* tmp = \frac{h^{1}}{1!}*A^{0} = I
       * and e_b = \Sum_{j=0}^{0}\frac{h^(j+1)}{(j+1)!}*A^{j} */
      for (size_type i=1; i< n; ++i) {
        /* tmp = \frac{h^{i}}{i!}*A^{i-1}
         * and e_b = \Sum_{j=0}^{i-1}\frac{h^(j+1)}{(j+1)!}*A^{j} */
        
        tmp *= (h/(i+1));
        tmp = prod(tmp,A);        
        /* tmp =  (h^(i+1)/(i+1)!)*A^i */

        e_b += tmp;
        /*  e_b = \Sum_{j=0}^{i}\frac{h^(j+1)}{(j+1)!}*A^{j} */
      }
      
      return (prod(e_b,b));
      /* out = ( \Sum_{j=0}^{n}\frac{h^(j+1)}{(j+1)!}*A^{j} ) b */
    }
    
    
    /* PAY ATTENTION!!! 
     * I supose that matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename Real>
    matrix<Real> 
    lu_decompose(const matrix<Real> &A, 
                 vector<size_type> &p_vect) 
    {
      Real max,sum,p_val;
      size_type i,j,k, size=A.size1(), pivot=0;
      
      vector<Real> scale(size);
      matrix<Real> O=A;
      
      for (i=0; i<size; i++) {
        max=0.0;
        for (j=0; j<size; j++) {
          if (abs(O(i,j)) > max ) {
            max=abs(O(i,j));
          }
        }
        if (max==0) {
          throw std::invalid_argument("The input matrix is singular");
        }
        scale(i)=1.0/max;
      }
      
      for (j=0; j<size; j++) {
        for (i=0; i<j; i++) {
          sum=O(i,j);
          for (k=0; k<i; k++) {
            sum -= O( i , k ) * O( k , j );
          }
          O(i,j)=sum;
        }
        
        max=0.0;
        
        for (i=j; i<size; i++) {
          sum=O(i,j);
          for (k=0; k<j; k++) {
            sum -= O( i , k ) * O( k , j );
          }
          O(i,j)=sum;
          p_val=sum*scale(i);
          if (abs(p_val) >= max) {
            max=abs(p_val);
            pivot=i;
          }
        }
        
        if (j != pivot ) {
          for (k=0; k<size; k++) {
            p_val=O(pivot,k);
            O(pivot,k)=O(j,k);
            O(j,k)=p_val;
          }
          scale(pivot)=scale(j);
        }
        
        p_vect(j)=pivot;
        
        if ( O(j,j) == 0 ) {
          throw std::invalid_argument("The input matrix is singular");
        }
    
        if (j < size-1) {
          for (i=j+1;i<size; i++) {
            O(i,j)=O(i,j)/O(j,j);
          }
        }
      }
      
      return O;
    }
    
    /* PAY ATTENTION!!! 
     * I supose that boost::numeric::ublas::matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename Real>
    vector<Real> 
    lu_solve(const matrix<Real> &A, 
             const vector<size_type> &p_vect, 
             const vector<Real> &b) 
    {
      size_type i_diag=0, idx_p, size=A.size1(),i,j;
      Real sum;
      vector<Real> sol=b;
    
      for (i=0; i<size; i++) {
        idx_p=p_vect(i);
        
        sum=sol(idx_p);
        sol(idx_p)=b(i);
        
        if (i_diag!=0) {
          for (j=i_diag-1; j<i; j++){
            
            sum -= (A( i , j ) * sol (j));
          }
        } else {
          if (sum!=0) {
            i_diag=i+1;
          }
        }
        
        sol(i)=sum;
        
      }
      
      for (idx_p=0; idx_p<size; idx_p++){
        i=(size-1)-idx_p;
        
        sum=sol(i);
        for (j=i+1; j< size; j++){
          sum -= A(i,j) * sol(j);
        }
        
        sol(i) = sum / A(i,i);
      }
      
      return sol;
    
    }
    
    template <typename Real>
    matrix<Real> 
    invert_matrix(const matrix<Real> &A) {
      
      size_type size=A.size1(),i,j;
      
      matrix<Real> inv_A(size,size);
      vector<Real> Id_vect(size), Id_sol;
      vector<size_type> p_vect(size);
      
      
      matrix<Real> lu_A=lu_decompose(A, p_vect);
      
      for (j=0; j<size; j++) {
        for (i=0; i<size; i++)
          Id_vect(i)=0.0;
        
        Id_vect(j)=1.0;
        
        Id_sol=lu_solve(lu_A,p_vect,Id_vect);
        
        for (i=0; i<size; i++)
          inv_A(i,j)=Id_sol(i);
      }
      
      return inv_A;
      
    }
    
    template <typename Real>
    inline Real find_matrix_denumerator(const matrix<Real> &A)
    {
      size_type i_max=A.size1(), j_max=A.size2();
      size_type i,j;
      
      if ((i_max==0)||(j_max==0)) 
        return 0.0;
      
      Real denum=denumerator(A(0,0));
      
      for (j=0; j< j_max; j++) {
        
        for (i=0; i< i_max; i++) {
          denum=lcm(numerator(denum),denumerator(A(i,j)));
        }
        
      }
      return denum;
    }
    
    
    template <typename Real>
    inline Real find_vector_denumerator(const vector<Real> &b) {
      
      size_type i_max=b.size();
      size_type i;
      
      if (i_max==0) 
        return 0.0;
      
      Real denum=denumerator(b(0));
      
      for (i=0; i< i_max; i++) {
        denum=lcm(numerator(denum),denumerator(b(i)));
      }
      
      return denum;
    }

  }
}

#endif /* _ARIADNE_LINEAR_ALGEBRA_H */
