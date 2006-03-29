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

namespace boost {
  namespace numeric {
    namespace ublas {

      template <typename Real>
      std::ostream&
      operator<<(std::ostream& os, const matrix<Real>& A)
      {
        if(A.size1()==0 || A.size2()==0) {
          return os << "[ ]";
        }
        
        for(uint i=0; i!=A.size1(); ++i) {
          os << (i==0 ? "[ " : "; ") << A(i,0);
          for(uint j=1; j!=A.size2(); ++j) {
            os << "," << A(i,j);
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

    template <typename Real>
    matrix<Real> zero_matrix(size_type r, size_type c) 
    {
      matrix<Real> A(r,c);
      for (size_type i=0; i<r; ++i) {
        for (size_type j=0; j<c; ++j) {
          A(i,j)=0.0;
        }
      }
      return A;
    }
  
    template<typename Real>
    Real
    norm(const matrix<Real>& A) 
    {
      Real result=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        Real row_sum=0;
        for(size_type j=0; j!=A.size2(); ++j) {
          row_sum+=abs(A(i,j));
        }
        result=std::max(result,row_sum);
      }
      return result;
    }
        
    
    
    template <typename Real>
    matrix<Real>
    exp_Ah(const matrix<Real>& A, 
           const Real& h, 
           const Real& e) 
    {
      matrix<Real> result=identity_matrix<Real>(A.size1());
      
      Real norm_Ah=h*norm(A);
      matrix<Real> AhpowNdivfN=result;
      uint n=0;
      while(norm(AhpowNdivfN)*n >= e*(n-norm_Ah)) {
        ++n;
        AhpowNdivfN=(h/n)*(AhpowNdivfN*A);
        result=result+AhpowNdivfN;
      }
      return result;
    }
    
    template <typename Real> 
    vector<Real> 
    exp_b_approx(const matrix<Real> &A, 
                 const vector<Real> &b, 
                 const Real h, 
                 const unsigned int n) 
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
   
    template <typename Real>
    void 
    lu_local_dec(matrix<Real> &A, 
                 const array<size_type> &row, const array<size_type> &col, 
                 const size_type &rows, const size_type &columns, 
                 const size_type &p) 
    {
      size_type i,j;
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
    matrix<Real> 
    lu_decompose(const matrix<Real> &A, 
                 array<size_type> &p_col, 
                 array<size_type> &p_row) 
    {
      //create output matrix
      matrix<Real> lu_A(A);
  
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
    template <typename Real>
    matrix<Real> 
    lu_decompose(const matrix<Real> &A, 
                 array<size_type> &p_array) 
    {
      Real max,sum,p_val;
      size_type i,j,k, rows=A.size1(), cols=A.size2(), pivot=0;
      
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
    template <typename Real>
    vector<Real> 
    lu_solve(const matrix<Real> &A, 
             const array<size_type> &p_array, 
             const vector<Real> &b) 
    {
      size_type i_diag=0, idx_p, rows=A.size1(), cols=A.size2(),i,j;
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

      size_type dim1=A.size1(),dim2=A.size2();
      vector<Real> x(dim1);
      matrix<Real> Q=identity_matrix<Real>(dim1),QA=A,Qi;
      Real norm,coef;

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
        
        Qi=identity_matrix<Real>(dim1);
          
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
   
    template <typename Real>
    matrix<Real>
    hermitian(const matrix<Real>& m) {
       return herm(m);
    }
    
    template <typename Real>
    matrix<Real> 
    inverse(const matrix<Real> &A) {
      
      size_type rows=A.size1(), cols=A.size2(),i,j;
      
      matrix<Real> inv_A(cols,rows);
      vector<Real> Id_vect(rows), Id_sol;
      array<size_type> p_array(rows);
      

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
    Integer 
    common_denominator(const matrix<Real>& A)
    {
      Integer denom=1;
      for (size_type i=0; i<A.size1(); ++i) {
        for (size_type j=0; j<A.size2(); ++j) {
          denom=lcm( denom, denominator(A(i,j)) );
        }
      }
      return denom;
    }
    
    template <typename Real>
    vector<Integer> 
    row_common_denominators(const matrix<Real>& A) 
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
    template <typename Real>
    void 
    transform_linear_inequalities(const matrix<Real>& T, 
                                  matrix<Real>& A, 
                                  vector<Real>& b) 
    {
      A=A*inverse(T);
    }
    
    template <>
    void 
    transform_linear_inequalities<Dyadic>(const matrix<Dyadic>& T, 
                                          matrix<Dyadic>& A, 
                                          vector<Dyadic>& b) 
    {
      typedef Dyadic Real;
      
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
      
      Real rmultiplier = convert_to<Real>(multiplier);
      matrix<Real> rTinv(n,n);
       for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          rTinv(i,j) = convert_to<Real>(iTinv(i,j));
        }
      }

      A=prod(A,rTinv);
      b=b*rmultiplier;
    }
   
    template <class Real>
    bool 
    independent_rows(matrix<Real> A) 
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
 
    template <typename Real>
    bool 
    have_same_dimensions(const matrix<Real> &A,  const matrix<Real> &B) 
    {
      return (A.size1()==B.size1() && A.size2()==B.size2());
    }
    
    template <typename Real>
    bool 
    equivalent_columns(const matrix<Real> &A, const size_type &A_col, 
                       const matrix<Real> &B, const size_type &B_col) 
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
  
    template<typename Real>
    size_type 
    find_first_not_null_in_col(const matrix<Real> &A, 
                               const size_type &col) 
    {
      size_type i=0;

      while ( (i<A.size1()) && (A(i,col)==0.0) ) {
        ++i;
      }
      return i;
      
    }
    
    template <class Real>
    matrix<Real> 
    remove_null_columns_but_one(const matrix<Real> &A) 
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

      matrix<Real> new_A(rows,std::max(1u,directions));

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
    
    template <typename Real>
    void 
    remove_null_columns(const matrix<Real> &A, 
                        array<size_type> &row, array<size_type> &col) 
    {
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


    template <typename Real>
    matrix<Real> 
    compute_space(const matrix<Real> &SA, 
                  array<size_type> &row,const array<size_type> &col) 
    {
      size_type cols=col.size(), rows=row.size();
   
      assert(cols>=rows);
   
      size_type SA_cols=SA.size2(), A_rows=SA_cols-rows;
      size_type i,j,k,j2;
   
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
