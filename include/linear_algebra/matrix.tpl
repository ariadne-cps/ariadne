/***************************************************************************
 *            matrix.tpl
 *
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
 
#include "../linear_algebra/matrix.h"

#include <algorithm>
#include "../base/array.h"

namespace boost {
  namespace numeric {
    namespace ublas {

      template <typename R>
      std::ostream&
      operator<<(std::ostream& os, const matrix<R>& A)
      {
        if(A.size1()==0 || A.size2()==0) {
          return os << "[ ]";
        }
        
        for(uint i=0; i!=A.size1(); ++i) {
          for(uint j=0; j!=A.size2(); ++j) {
            os << (j==0 ? (i==0 ? "[ " : "; ") : ",");
            os << A(i,j);
          }
        }
        os << " ]";
        return os;
      }
       
      template <>
      std::ostream&
      operator<<(std::ostream& os, const matrix<Ariadne::Rational>& A)
      {
        for(uint i=0; i!=A.size1(); ++i) {
          for(uint j=0; j!=A.size2(); ++j) {
            os << (j==0 ? (i==0 ? "[ " : "; ") : ",");
            os << Ariadne::convert_to<double>(A(i,j));
          }
        }
        os << " ]";
        return os;
      }
       

    }
  }
}
      

namespace Ariadne {
  namespace LinearAlgebra {

    template <typename R>
    matrix<R> zero_matrix(size_type r, size_type c) 
    {
      matrix<R> A(r,c);
      for (size_type i=0; i<r; ++i) {
        for (size_type j=0; j<c; ++j) {
          A(i,j)=0.0;
        }
      }
      return A;
    }
  
    template<typename R>
    R
    norm(const matrix<R>& A) 
    {
      R result=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        R row_sum=0;
        for(size_type j=0; j!=A.size2(); ++j) {
          row_sum+=abs(A(i,j));
        }
        result=std::max(result,row_sum);
      }
      return result;
    }
        
    template<typename R>
    R
    log_norm(const matrix<R>& A) 
    {
      R result=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        R row_sum=A(i,i);
        for(size_type j=0; j!=A.size2(); ++j) {
          if(i!=j) {
            row_sum+=abs(A(i,j));
          }
        }
        result=max(result,row_sum);
      }
      return result;
    }
        
    template<typename R>
    matrix<R>
    concatenate_columns(const matrix<R>& A1, const matrix<R>& A2) {
      assert(A1.size1()==A2.size1());
      LinearAlgebra::matrix<R> result(A1.size1(),A1.size2()+A2.size2());
      for(size_type i=0; i!=result.size1(); ++i) {
        for(size_type j=0; j!=A1.size2(); ++j) {
          result(i,j)=A1(i,j);
        }
        for(size_type j=0; j!=A2.size2(); ++j) {
          result(i,j+A1.size2())=A2(i,j);
        }
      }
      return result;
    }
    

   
    
    template <typename R>
    matrix<R>
    exp_approx(const matrix<R>& A, 
               const R& e) 
    {
      matrix<R> result=identity_matrix<R>(A.size1());
      
      R norm_A=norm(A);
      matrix<R> term=result;
      uint n=0;
      while(norm(term)*n >= e*(n-norm_A)) {
        ++n;
        term=(term*A)/n;
        result=result+term;
      }
      return result;
    }
    

    template <typename R>
    void 
    lu_local_dec(matrix<R> &A, 
                 const array<size_type> &row, const array<size_type> &col, 
                 const size_type &rows, const size_type &columns, 
                 const size_type &p) 
    {
      size_type i,j;
      R coef;
  
      // perform lu decomposition on the sub matrix
      for (i=p+1; i< rows; i++) {
        coef=A(row[i],col[p])/A(row[p],col[p]);
    
        for (j=p+1; j< columns; j++) 
          A(row[i],col[j])-=(A(row[p],col[j])*coef);
    
        A(row[i],col[p])=coef;
    
      }
    }
   
    template <typename R>
    matrix<R> 
    lu_decompose(const matrix<R> &A, 
                 array<size_type> &p_col, 
                 array<size_type> &p_row) 
    {
      using std::swap;
      
      //create output matrix
      matrix<R> lu_A(A);
  
      size_type rows,columns;
      rows=A.size1();
      columns=A.size2();

      // create the two pivot vectors
      array<size_type> row(rows);
      p_col=array<size_type>(columns);

      size_type i,j;
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
 
      p_row=array<size_type>(rows);
  
      for (i=0; i<rows; i++) 
        p_row[i]=row[i];
    
      return lu_A;
    }
    
    /* PAY ATTENTION!!! 
     * I supose that matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename R>
    matrix<R> 
    lu_decompose(const matrix<R> &A, 
                 array<size_type> &p_array) 
    {
      R max,sum,p_val;
      size_type i,j,k, rows=A.size1(), cols=A.size2(), pivot=0;
      
      vector<R> scale(rows);
      matrix<R> O=A;
      
      for (i=0; i<rows; i++) {
        max=0.0;
        for (j=0; j<cols; j++) {
          if (abs(O(i,j)) > max ) {
            max=abs(O(i,j));
          }
        }
        if (max==0) {
          throw std::runtime_error("Matrix is singular");
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
          throw std::runtime_error("Matrix is singular");
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
    template <typename R>
    vector<R> 
    lu_solve(const matrix<R> &A, 
             const array<size_type> &p_array, 
             const vector<R> &b) 
    {
      size_type i_diag=0, idx_p, rows=A.size1(), cols=A.size2(),i,j;
      R sum;
      vector<R> sol(cols);
   
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
    template <typename R>
    inline
    matrix<R> 
    Householder_QR(const matrix<R> &A) {

      size_type dim1=A.size1(),dim2=A.size2();
      vector<R> x(dim1);
      matrix<R> Q=identity_matrix<R>(dim1),QA=A,Qi;
      R norm,coef;

      for (size_type i=0; i< dim2; i++) {
        
        norm=0;
        
        for (size_type j=i; j< dim1; j++) {
          norm+=QA(j,i)*QA(j,i);
          x(j)=QA(j,i);
        }
        coef=norm-x(i)*x(i);
        
        if (x(i)>=0)  
           x(i)-=sqrt(norm);
        else 
           x(i)+=sqrt(norm);
      
        coef+=x(i)*x(i);
        
        Qi=identity_matrix<R>(dim1);
          
        if (coef!=0) {
          for (size_type j=i; j< dim2; j++) {
            for (size_type k=i; k< dim1; k++) {
               Qi(k,j)-=((2*x(k)*x(j))/coef);
            }
          }
        }
        
        QA=prod(Qi,QA);
        Q=prod(Q,Qi);
      }
     
      return Q;
    }
   
    template <typename R>
    matrix<R>
    hermitian(const matrix<R>& m) {
       return herm(m);
    }
    
    template <typename R>
    matrix<R> 
    inverse(const matrix<R> &A) {
      
      size_type rows=A.size1(), cols=A.size2(),i,j;
      
      matrix<R> inv_A(cols,rows);
      vector<R> Id_vect(rows), Id_sol;
      array<size_type> p_array(rows);
      

      matrix<R> lu_A=lu_decompose(A, p_array);
     
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
    
    template <typename R>
    Integer 
    common_denominator(const matrix<R>& A)
    {
      Integer denom=1;
      for (size_type i=0; i<A.size1(); ++i) {
        for (size_type j=0; j<A.size2(); ++j) {
          denom=lcm( denom, denominator(A(i,j)) );
        }
      }
      return denom;
    }
    
    template <typename R>
    vector<Integer> 
    row_common_denominators(const matrix<R>& A) 
    {
      vector<Integer> denoms(A.size1());
      for(size_type i=0; i!=A.size1(); ++i) {
        Integer denom=1;
        for(size_type j=0; j!=A.size2(); ++j) {
          denom=lcm( denom, denominator(A(i,j)) );
        }
        denoms(i)=denom;
      }
      return denoms;
    }

    
    
    /* \brief Transforms the linear inequalities $Ax\leq b$ to $AT^{-1}y \leq b$. */
    template <typename R>
    void 
    transform_linear_inequalities(const matrix<R>& T, 
                                  matrix<R>& A, 
                                  vector<R>& b) 
    {
      A=A*inverse(T);
    }
    
    template <>
    void 
    transform_linear_inequalities<Dyadic>(const matrix<Dyadic>& T, 
                                          matrix<Dyadic>& A, 
                                          vector<Dyadic>& b) 
    {
      typedef Dyadic R;
      
      size_type n=A.size1();
      if(A.size1() != b.size()) {
        throw std::domain_error("Invalid linear inequalities"); 
      }
      if(n!=T.size1() || n!=T.size2()) {
        throw std::domain_error("Invalid linear transformation");
      }
      
      matrix<Rational> Trat(T.size1(),T.size2());
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Trat(i,j)=convert_to<Rational>(T(i,j));
        }
      }
      matrix<Rational> Tinv=inverse(Trat);
      Trat.clear();
      
      Integer multiplier=common_denominator(Tinv);
      
      matrix<Integer> iTinv(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          iTinv(i,j) = numerator(Tinv(i,j)) * (multiplier/denominator(Tinv(i,j)));
        }
      }
      
      R rmultiplier = convert_to<R>(multiplier);
      matrix<R> rTinv(n,n);
       for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          rTinv(i,j) = convert_to<R>(iTinv(i,j));
        }
      }

      A=prod(A,rTinv);
      b=b*rmultiplier;
    }
   
    template <class R>
    bool 
    independent_rows(matrix<R> A) 
    {
     const size_type rows = A.size1(),cols = A.size2();
     size_type i,j,i2,j2;

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
 
    template <typename R>
    bool 
    have_same_dimensions(const matrix<R> &A,  const matrix<R> &B) 
    {
      return (A.size1()==B.size1() && A.size2()==B.size2());
    }
    
    template <typename R>
    bool 
    equivalent_columns(const matrix<R> &A, const size_type &A_col, 
                       const matrix<R> &B, const size_type &B_col) 
    {
      if (A.size1()!=B.size1()) {
        throw std::domain_error("The two matrix have a diffentent number of rows"); 
      }
      for (size_type i=0; i<A.size1(); i++) {
        if (A(i,A_col)!=B(i,B_col)) {
          return false;
        }
      }
      return true;
    }
  
    template<typename R>
    size_type 
    find_first_not_null_in_col(const matrix<R> &A, 
                               const size_type &col) 
    {
      size_type i=0;

      while ( (i<A.size1()) && (A(i,col)==0.0) ) {
        ++i;
      }
      return i;
      
    }
    
    template <class R>
    matrix<R> 
    remove_null_columns_but_one(const matrix<R> &A) 
    {
      size_type directions=0;
      size_type rows=A.size1();
      size_type cols=A.size2();
      array<bool> not_null(cols);

      if (cols<2) {
        return A;
      }

      for(size_type j=0; j!=cols; ++j) {
        not_null[j]=(find_first_not_null_in_col(A,j)<rows);
        if (not_null[j]) {
          directions++; 
        }
      }

      matrix<R> new_A(rows,std::max(1u,directions));

      size_type j2=0;
      for(size_type j=0; j!=cols; ++j) {
        if (not_null[j]) { 
          for(size_type i=0; i!=rows; ++i) {
            new_A(i,j2)=A(i,j);
          }
          j2++;
        }
      }

      return new_A;
    }
    
    template <typename R>
    void 
    remove_null_columns(const matrix<R> &A, 
                        array<size_type> &row, array<size_type> &col) 
    {
      using std::swap;
      
      size_type i,j, columns=col.size(), rows=row.size();
      size_type new_columns=columns;
  
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
        array<size_type> p_col(new_columns);

        for (size_type j=0; j<new_columns; j++) 
          p_col[j]=col[j];

        col=p_col;
      }
    }


    template <typename R>
    matrix<R> 
    compute_space(const matrix<R> &SA, 
                  array<size_type> &row,const array<size_type> &col) 
    {
      size_type cols=col.size(), rows=row.size();
   
      assert(cols>=rows);
   
      size_type SA_cols=SA.size2(), A_rows=SA_cols-rows;
      size_type i,j,k,j2;
   
      matrix<R> A(A_rows,SA_cols);

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

      R aux;
  
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
