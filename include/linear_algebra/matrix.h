/***************************************************************************
 *            matrix.h
 *
 *  Mon May  3 12:31:15 2004
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
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
 
/*! \file matrix.h
 *  \brief Matrices and matrix operations.
 */

#ifndef _ARIADNE_MATRIX_H
#define _ARIADNE_MATRIX_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include "../base/basic_type.h"
#include "../base/numerical_type.h"
#include "../base/interval.h"

namespace Ariadne {
  namespace LinearAlgebra {

    using boost::numeric::ublas::identity_matrix;
    using boost::numeric::ublas::herm;
    using boost::numeric::ublas::vector;
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::matrix_column;
  
    class SingularMatrix{
    public:
       SingularMatrix(){}
       ~SingularMatrix(){}    
    };
    
    template <typename Real>
    matrix< Interval<Real> > 
    prod(const matrix< Interval<Real> >& A, const matrix<Real>& B) {
      matrix< Interval<Real> > result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=Interval<Real>(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }
      
    template <typename Real>
    matrix< Interval<Real> > 
    prod(const matrix<Real>& A, const matrix< Interval<Real> >& B) {
      matrix< Interval<Real> > result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=Interval<Real>(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }

    template <typename Real>
    inline matrix<Real> zero_matrix(dimension_type dim) {
      matrix<Real> A(dim,dim);
      for (dimension_type i=0; i<dim; ++i) {
        for (dimension_type j=0; j<dim; ++j) {
          A(i,j)=0.0;
        }
      }
      return A;
    }
  
    /*
    template <typename Real>
    inline matrix<Real> identity_matrix(dimension_type dim) {
      matrix<Real> A(dim,dim);
      for (dimension_type i=0; i<dim; i++) {
        for (dimension_type j=0; j<dim; j++) {
          A(i,j)=0.0;
        }
        A(i,i)=1.0;
      }
      return A;
    }
    */
    
    template <typename Real>
    inline matrix<Real> exp_Ah(const matrix<Real> &A, 
                               const Real h, const unsigned int n) 
    {
      matrix<Real> tmp,e_Ah;
      
      e_Ah=identity_matrix<Real>(A.size1());
      tmp=e_Ah;
      
      /* tmp = \frac{h^{0}}{0!}*A^{0} = I
       * and e_Ah = \Sum_{j=0}^{0}\frac{h^j}{j!}*A^{j} = I */
      for (dimension_type i=1; i<n; ++i) {
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
      for (dimension_type i=1; i< n; ++i) {
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
                 vector<dimension_type> &p_vect) 
    {
      Real max,sum,p_val;
      dimension_type i,j,k, rows=A.size1(), cols=A.size2(), pivot=0;
      
      vector<Real> scale(rows);
      matrix<Real> O=A;
      
      for (i=0; i<rows; i++) {
        max=0.0;
        for (j=0; j<cols; j++) {
          if (abs(O(i,j)) > max ) {
            max=abs(O(i,j));
          }
        }
        if (max==0) {
          throw SingularMatrix();
        }
        scale(i)=1.0/max;
      }
      
      for (j=0; j<cols && j<rows; j++) {
        for (i=0; i<j; i++) {
          sum=O(i,j);
          for (k=0; k<i; k++) {
            sum -= O( i , k ) * O( k , j );
          }
          O(i,j)=sum;
        }
        
        max=0.0;
        
        for (i=j; i<rows; i++) {
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
          for (k=0; k<cols; k++) {
            p_val=O(pivot,k);
            O(pivot,k)=O(j,k);
            O(j,k)=p_val;
          }
          scale(pivot)=scale(j);
        }
      
        p_vect(j)=pivot;
        
        if ( O(j,j) == 0 ) {
          throw SingularMatrix();
        }
    
        if (j < cols-1) {
          for (i=j+1;i<rows; i++) {
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
             const vector<dimension_type> &p_vect, 
             const vector<Real> &b) 
    {
      dimension_type i_diag=0, idx_p, rows=A.size1(), cols=A.size2(),i,j;
      Real sum;
      vector<Real> sol(cols);
   
      for (i=0; i<rows&& i<cols; i++) {
        sol(i)=b(i);
      }
      
      for (i=0; i<rows; i++) {
        idx_p=p_vect(i);
       
	if (idx_p < cols) {
          sum=sol(idx_p);
          sol(idx_p)=b(i);
        
          if (i_diag!=0) {
            for (j=i_diag-1; j<i; j++){
              if (j< cols) sum -= (A( i , j ) * sol (j));
            }
          } else {
            if (sum!=0) {
              i_diag=i+1;
            }
          }
        
          if (i< cols) sol(i)=sum;
        }  
      }
      
      for (idx_p=0; idx_p<rows; idx_p++){
        i=(rows-1)-idx_p;
        
        if ((i<cols)&&(i<rows)) {
	  sum=sol(i);
          for (j=i+1; j< cols; j++){
            sum -= A(i,j) * sol(j);
          }
        
          sol(i) = sum / A(i,i);
	}
      }
      
      return sol;
    
    }
  
    /* WARNING!!! The following function has some precision problems */ 
    template <typename Real>
    inline
    matrix<Real> 
    Householder_QR(const matrix<Real> &A) {

      dimension_type dim1=A.size1(),dim2=A.size2();
      vector<Real> x(dim1);
      matrix<Real> Q=identity_matrix<Real>(dim1),QA=A,Qi;
      Real norm,coef;

      for (dimension_type i=0; i< dim2; i++) {
	
	norm=0;
	
	for (dimension_type j=i; j< dim1; j++) {
	  norm+=QA(j,i)*QA(j,i);
	  x(j)=QA(j,i);
	}
	coef=norm-x(i)*x(i);
	
	if (x(i)>=0)  
	   x(i)-=sqrt(norm);
	else 
	   x(i)+=sqrt(norm);

	coef+=x(i)*x(i);
	
	Qi=identity_matrix<Real>(dim1);
	  
	if (coef!=0) {
	  for (dimension_type j=i; j< dim2; j++) {
	    for (dimension_type k=i; k< dim1; k++) {
	       Qi(k,j)-=((2*x(k)*x(j))/coef);
	    }
	  }
	}
	
	QA=prod(Qi,QA);
	Q=prod(Q,Qi);
      }
     
      return Q;
    }
   
    template <typename Real>
    inline
    matrix<Real>
    hermitian(const matrix<Real>& m) {
       return herm(m);
    }
    
    template <typename Real>
    inline
    dimension_type
    number_of_rows(const matrix<Real> &A) {
      return A.size1();
    }
    
    template <typename Real>
    inline
    dimension_type
    number_of_columns(const matrix<Real> &A) {
      return A.size2();
    }
    
    template <typename Real>
    inline
    matrix<Real> 
    inverse(const matrix<Real> &A) {
      
      dimension_type rows=A.size1(), cols=A.size2(),i,j;
      
      matrix<Real> inv_A(cols,rows);
      vector<Real> Id_vect(rows), Id_sol;
      vector<dimension_type> p_vect(rows);
      

      matrix<Real> lu_A=lu_decompose(A, p_vect);
     
      for (j=0; j<rows; j++) {
        for (i=0; i<rows; i++)
          Id_vect(i)=0.0;
        
        Id_vect(j)=1.0;

        Id_sol=lu_solve(lu_A,p_vect,Id_vect);
        
        for (i=0; i<cols; i++)
          inv_A(i,j)=Id_sol(i);
      }
      
      return inv_A;
    }
    
    /* FIXME:Disable this template since the inverse of a dyadic matrix cannot be computed. */
    /*
    template < >
    matrix<Dyadic>
    inverse(const matrix<Dyadic>& A);
    */
    
    template <typename Real>
    matrix<Real> 
    invert_matrix(const matrix<Real> &A) {
      return inverse(A);
    }
    
    template <typename Real>
    inline Integer common_denominator(const matrix<Real>& A)
    {
      Integer denom=1;
      for (dimension_type i=0; i<A.size1(); ++i) {
        for (dimension_type j=0; j<A.size2(); ++j) {
          denom=lcm( denom, denominator(A(i,j)) );
        }
      }
      return denom;
    }
    
    template <typename Real>
    inline vector<Integer> row_common_denominators(const matrix<Real>& A) 
    {
      vector<Integer> denoms(A.size1());
      for(dimension_type i=0; i!=A.size1(); ++i) {
        Integer denom=1;
        for(dimension_type j=0; j!=A.size2(); ++j) {
          denom=lcm( denom, denominator(A(i,j)) );
        }
        denoms(i)=denom;
      }
      return denoms;
    }

    
    
    /* \brief Transforms the linear inequalities $Ax\leq b$ to $At^{-1}y \leq b$. */
    template <typename R>
    inline void transform_linear_inequalities(const matrix<R>& T, 
                                       matrix<R>& A, 
                                       vector<R>& b) 
    {
      A=A*inverse(T);
    }
    
    template <>
    inline void transform_linear_inequalities<Dyadic>(const matrix<Dyadic>& T, 
                                               matrix<Dyadic>& A, 
                                               vector<Dyadic>& b) 
    {
      typedef Dyadic Real;
      
      dimension_type n=A.size1();
      if(A.size1() != b.size()) {
        throw std::domain_error("Invalid linear inequalities"); 
      }
      if(n!=T.size1() || n!=T.size2()) {
        throw std::domain_error("Invalid linear transformation");
      }
      
      matrix<Rational> Trat(T.size1(),T.size2());
      for(dimension_type i=0; i!=n; ++i) {
        for(dimension_type j=0; j!=n; ++j) {
          Trat(i,j)=convert_to<Rational>(T(i,j));
        }
      }
      matrix<Rational> Tinv=inverse(Trat);
      Trat.clear();
      
      Integer multiplier=common_denominator(Tinv);
      
      matrix<Integer> iTinv(n,n);
      for(dimension_type i=0; i!=n; ++i) {
        for(dimension_type j=0; j!=n; ++j) {
          iTinv(i,j) = numerator(Tinv(i,j)) * (multiplier/denominator(Tinv(i,j)));
        }
      }
      
      Real rmultiplier = convert_to<Real>(multiplier);
      matrix<Real> rTinv(n,n);
       for(dimension_type i=0; i!=n; ++i) {
        for(dimension_type j=0; j!=n; ++j) {
          rTinv(i,j) = convert_to<Real>(iTinv(i,j));
        }
      }

      A=prod(A,rTinv);
      b=b*rmultiplier;
    }
   
    template <class T>
    inline bool independent_rows(matrix<T> A) {
     const size_t rows = A.size1(),cols = A.size2();
     size_t i,j,i2,j2;

      if (rows > cols) return false;
    
      for ( i=0; i<rows; i++) {
        j=0;
        while (A(i,j)==0) {
          j++; 
          if (j >= cols) return false;
        }
        for (i2=i+1; i2<rows; i2++) {
          for (j2=0; j2<cols; j2++) {
            A(i2,j2)= A(i2,j2) - ((A(i,j2) * A(i2,j))/A(i,j));
          }
        }
      }

      return true;
    } 
   
    template <class T>
    inline
    matrix<T> remove_null_columns_but_one(const matrix<T> &A) {
      size_t point=0,directions=0;
      size_t cols=number_of_columns(A), 
      rows=number_of_rows(A);
      array<bool> null_v(cols);

      if (cols==0) {
        return A;
      }

      for(size_t j=0; j!=cols; ++j) {
        null_v[j]=false;
        for(size_t i=0; i!=rows; ++i) {
          if (A(i,j)!=0.0) { null_v[j]=true; }
        }
        if (null_v[j]) { directions++; }
      }

      matrix<T> new_A(rows,std::max((size_t)1,directions));

      for(size_t j=0; j!=cols; ++j) {
        if (null_v[j]) { 
          for(size_t i=0; i!=rows; ++i) {
            new_A(i,point)=A(i,j);
          }
          ++point;
        }
      }

      return new_A;
    }

  }
}

namespace Ariadne {
  namespace Evaluation {

    //FIXME: Hack to include matrix/vector product operators.
    template<typename R>
    inline
    LinearAlgebra::vector<R>
    operator*(const LinearAlgebra::matrix<R>& A, const LinearAlgebra::vector<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    LinearAlgebra::matrix<R>
    operator*(const LinearAlgebra::matrix<R>& A, const LinearAlgebra::matrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    LinearAlgebra::matrix< Interval<R> >
    operator*(const LinearAlgebra::matrix< Interval<R> >& A, const LinearAlgebra::matrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    LinearAlgebra::matrix< Interval<R> >
    operator*(const LinearAlgebra::matrix<R>& A, const LinearAlgebra::matrix< Interval<R> >& B) {
      return prod(A,B);
    }
 

  }
}

namespace Ariadne {
  namespace Geometry {

    //FIXME: Hack to include matrix/vector product operators.
    template<typename R>
    inline
    LinearAlgebra::vector<R>
    operator*(const LinearAlgebra::matrix<R>& A, const LinearAlgebra::vector<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    LinearAlgebra::matrix<R>
    operator*(const LinearAlgebra::matrix<R>& A, const LinearAlgebra::matrix<R>& B) {
      return prod(A,B);
    }
  }
}


namespace boost { namespace numeric { namespace ublas {

    template <typename Real>
    std::ostream&
    operator<<(std::ostream& os, const matrix<Real>& A)
    {
      if(A.size1()==0 || A.size2()==0) {
        return os << "[ ]";
      }
      
      for(uint i=0; i!=A.size1(); ++i) {
        os << (i==0 ? "[ " : " ; ") << A(i,0);
        for(uint j=1; j!=A.size2(); ++j) {
          os << "," << A(i,j);
        }
      }
      os << " ]";
      return os;
    }
    
}}}


#endif /* _ARIADNE_MATRIX_H */
