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
#include "../base/utility.h"

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
   
   template <typename Real>
   inline 
   void lu_local_dec(matrix<Real> &A, const array<size_t> &row, 
		const array<size_t> &col, const size_t &rows, 
		const size_t &columns, const size_t &p) {
      size_t i,j;
      Real coef;
	
      // perform lu decomposition on the sub matrix
      for (i=p+1; i< rows; i++) {
        coef=A(row[i],col[p])/A(row[p],col[p]);
	  
        for (j=p+1; j< columns; j++) 
          A(row[i],col[j])-=(A(row[p],col[j])*coef);
    
        A(row[i],col[p])=coef;
    
      }
    }
   
    template <typename Real>
    inline	    
    matrix<Real> lu_decompose(const matrix<Real> &A, 
        array<size_t> &p_col, array<size_t> &p_row) {
	
      //create output matrix
      matrix<Real> lu_A(A);
	
      size_t rows,columns;
      rows=number_of_rows(A);
      columns=number_of_columns(A);

      // create the two pivot vectors
      array<size_t> row(rows);
      p_col=array<size_t>(columns);

      size_t i,j;
      for (j=0; j< columns; j++) 
        p_col[j]=j;

      for (i=0; i< rows; i++) 
        row[i]=i;
  
      i=0;
      // for all the linear independent rows
      while (i<rows) {
        j=i;

        //search a not null element in the row
        while ((j<columns)&&(lu_A(row[i],p_col[j])==0.0)) 
          j++;
    
        // if it does not exist, decrease the linear independent row number and 
        // swap the last linear independent row with the current one
        if (j==columns) {
          rows--;
          swap(row[i],row[rows]);
        } else { // otherwise (if the current row is linear independent for the 
	         // previous ones)
      
          // swap the current column with the first column which has a not null 
          // value in the i-th row
          swap(p_col[i],p_col[j]);

          // perform lu decomposition on the sub matrix A(i..rows,i..columns)
          lu_local_dec(lu_A, row, p_col, rows, columns, i);
       
          i++;
        }
      }
 
      p_row=array<size_t>(rows);
  
      for (i=0; i<rows; i++) 
        p_row[i]=row[i];
    
      return lu_A;
    }
    
    /* PAY ATTENTION!!! 
     * I supose that matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename Real>
    inline
    matrix<Real> 
    lu_decompose(const matrix<Real> &A, 
                 array<dimension_type> &p_array) 
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
      
        p_array[j]=pivot;
        
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
             const array<dimension_type> &p_array, 
             const vector<Real> &b) 
    {
      dimension_type i_diag=0, idx_p, rows=A.size1(), cols=A.size2(),i,j;
      Real sum;
      vector<Real> sol(cols);
   
      for (i=0; i<rows&& i<cols; i++) {
        sol(i)=b(i);
      }
      
      for (i=0; i<rows; i++) {
        idx_p=p_array[i];
       
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
      array<dimension_type> p_array(rows);
      

      matrix<Real> lu_A=lu_decompose(A, p_array);
     
      for (j=0; j<rows; j++) {
        for (i=0; i<rows; i++)
          Id_vect(i)=0.0;
        
        Id_vect(j)=1.0;

        Id_sol=lu_solve(lu_A,p_array,Id_vect);
        
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

    
    
    /* \brief Transforms the linear inequalities $Ax\leq b$ to $AT^{-1}y \leq b$. */
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
 
    template <typename Real>
    inline 
    bool have_same_dimensions(const matrix<Real> &A,  const matrix<Real> &B) {

      return ((number_of_columns(A)==number_of_columns(B))&&
	                (number_of_rows(A)==number_of_rows(B)));
    }
    
    template <typename Real>
    inline 
    bool equivalent_columns(const matrix<Real> &A, 
        const size_t &A_col, const matrix<Real> &B, 
        const size_t &B_col) {

      if (number_of_rows(A)!=number_of_rows(B))
        throw std::domain_error("The two matrix have a diffentent number of rows"); 

      for (size_t i=0; i< number_of_rows(A); i++)
	if (A(i,A_col)!=B(i,B_col)) return false;
   
      return true;
    }
  
    template<typename R>
    inline 
    size_t find_first_not_null_in_col(const matrix<R> &A, 
		    const size_t &col) {
	    
      size_t i=0;

      while ((i< number_of_rows(A))&&(A(i,col)==0.0))
      	i++;

      return i;
      
    }
    
    template <class T>
    inline
    matrix<T> remove_null_columns_but_one(const matrix<T> &A) {
      size_t directions=0;
      size_t cols=number_of_columns(A), rows=number_of_rows(A);
      array<bool> not_null(cols);

      if (cols<2) {
        return A;
      }

      for(size_t j=0; j!=cols; ++j) {
        not_null[j]=(find_first_not_null_in_col(A,j)<rows);
        if (not_null[j]) 
	  directions++; 
      }

      matrix<T> new_A(rows,std::max((size_t)1,directions));

      size_t j2=0;
      for(size_t j=0; j!=cols; ++j) {
        if (not_null[j]) { 
          for(size_t i=0; i!=rows; ++i) {
            new_A(i,j2)=A(i,j);
          }
          j2++;
        }
      }

      return new_A;
    }
    
    template <typename Real>
    inline 
    void remove_null_columns(const matrix<Real> &A, 
        array<size_t> &row, array<size_t> &col) {

      size_t i,j, columns=col.size(), rows=row.size();
      size_t new_columns=columns;
  
      j=columns-1;
      while (j>rows) {
        i=0;
        while ((i<rows)&&(A(row[i],col[j])==0)) 
          i++;
    
        if (i==rows) {
          swap(col[j],col[rows+1]);
          new_columns--;
        }
    
        j--;
      }
  
      if (new_columns<columns) {
        array<size_t> p_col(new_columns);

        for (size_t j=0; j<new_columns; j++) 
          p_col[j]=col[j];

        col=p_col;
      }
    }


    template <typename Real>
    inline 
    matrix<Real> compute_space(const matrix<Real> &SA, 
        array<size_t> &row,const array<size_t> &col) {
	
      size_t cols=col.size(), rows=row.size();
   
      assert(cols>=rows);
	   
      size_t SA_cols=number_of_columns(SA), A_rows=SA_cols-rows;
      size_t i,j,k,j2;
   
      matrix<Real> A(A_rows,SA_cols);

      k=0;
      for (i=0; i<SA_cols; i++) {
        j=0;
        while ((j<cols)&&(col[j]!=i))
          j++;

        if (j==cols) {
          A(k,i)=1;
          k++;
        }
      }

      Real aux;
  
      for (i=rows; i<cols; i++) {
        A(k,col[i])=1;
   
        A(k,col[rows-1])=-SA(row[rows-1], col[i])/SA(row[rows-1],col[rows-1]);

        j=rows; 
        while (j>0) {
          j--;
       
          aux=-SA(row[j], col[i]);
          for (j2=j+1; j2<rows; j2++)
            aux-=SA(row[j], col[j2])*A(k,col[j2]);
 
          A(k,col[j])=aux/SA(row[j], col[j]);
        }
      }

      return A;
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
