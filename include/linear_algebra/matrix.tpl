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

#include "matrix.h"

#include <algorithm>
#include "../base/array.h"

#include "../numeric/integer.h"
#include "../numeric/arithmetic.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/lu_matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    template <typename R>
    Matrix<R>::Matrix(const std::string& s)
    {
        std::stringstream ss(s);
        ss >> *this;
    }

    template<typename R>
    bool
    Matrix<R>::operator==(const Matrix<R>& A) const
    {
      if(this->number_of_rows()!=A.number_of_rows() ||
            this->number_of_columns()!=A.number_of_columns())
      {
        return false;
      }
      for(size_type i=0; i!=this->number_of_rows(); ++i) {
        for(size_type j=0; j!=this->number_of_columns(); ++j) {
          if((*this)(i,j)!=A(i,j)) { return false; }
        }
      }
      return true;
    }

    template<typename R>
    bool
    Matrix<R>::operator!=(const Matrix<R>& A) const
    {
      return !(*this==A);
    }

    template<typename R>
    R
    Matrix<R>::norm() const
    {
      const Matrix<R>& A=*this;
      F result=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        F row_sum=0;
        for(size_type j=0; j!=A.size2(); ++j) {
          row_sum+=abs(A(i,j));
        }
        result=max(result,row_sum);
      }
      return result;
    }

    template<typename R>
    R
    Matrix<R>::log_norm() const
    {
      const Matrix<R>& A=*this;
      F result=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        F row_sum=A(i,i);
        for(size_type j=0; j!=A.size2(); ++j) {
          if(i!=j) {
            row_sum+=abs(A(i,j));
          }
        }
        result=max(result,row_sum);
      }
      return result;
    }

    template <typename R>
    Matrix<R>
    Matrix<R>::transpose() const
    {
      size_type m=this->number_of_rows();
      size_type n=this->number_of_columns();
      const Matrix<R>& self=*this;

      Matrix<R> result(n,m);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=m; ++j) {
          result(i,j)=self(j,i);
        }
      }
      return result;
    }


    template <typename R>
    bool
    Matrix<R>::singular() const {
      return LUMatrix<R>(*this).singular();
    }

    template <typename R>
    R
    Matrix<R>::determinant() const
    {
      return LUMatrix<R>(*this).determinant();
    }

    template<typename R>
    Matrix<R>
    Matrix<R>::concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2) {
      assert(A1.size1()==A2.size1());
      LinearAlgebra::Matrix<R> result(A1.size1(),A1.size2()+A2.size2());
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
    Matrix<R>
    Matrix<R>::inverse() const
    {
      return LUMatrix<R>(*this).inverse();
    }


    template <typename R>
    Vector<R>
    Matrix<R>::solve(const Vector<R>& v) const
    {
      return this->inverse()*v;
    }


    template <typename R>
    std::ostream&
    Matrix<R>::write(std::ostream& os) const
    {
      const Matrix<R>& A=*this;
      os << "[";
      for(uint i=0; i!=A.size1(); ++i) {
        for(uint j=0; j!=A.size2(); ++j) {
          os << (j==0 ? (i==0 ? " " : "; ") : ",");
          //os << Ariadne::convert_to<double>(A(i,j));
          os << A(i,j);
        }
      }
      os << " ]";
      return os;
    }


    template <typename R>
    std::istream&
    Matrix<R>::read(std::istream& is)
    {
      char c;
      is >> c;
      is.putback(c);
      if(c=='[') {
        is >> c;
        /* Representation as a literal [a11,a12,...,a1n; a21,a22,...a2n; ... ; am1,am2,...,amn] */
        std::vector< std::vector<R> > v;
        R x;
        c=';';
        while(is && c==';') {
          v.push_back(std::vector<R>());
          c=',';
          while(is && c==',') {
            is >> x;
            v.back().push_back(x);
            is >> c;
          }
        }
        if(is) {
          Matrix<R>& A=*this;
          A=Matrix<R>(v.size(),v.front().size());
          for(size_type i=0; i!=A.size1(); ++i) {
            assert(v[i].size()==A.size2());
            for(size_type j=0; j!=A.size2(); ++j) {
              A(i,j)=v[i][j];
            }
          }
        }
      }
      else {
        assert(false);
      }
      return is;
    }





    // FIXME:
    // The following functions are only used by the Zonotope class.
    // Since they are not very general, we should probably remove 
    // many of these operations.
    // Also note that many algorithms can probably be better implemented
    // by using LU or QR factorizations.

    template <class R>
    bool
    independent_rows(Matrix<R> A)
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
            A(i2,j2)=(A(i2,j2)*A(i,j) - A(i,j2)*A(i2,j));
          }
        }
      }

      return true;
    }


    template <typename R>
    bool
    equivalent_columns(const Matrix<R> &A, const size_type &A_col,
                       const Matrix<R> &B, const size_type &B_col)
    {
      if (A.size1()!=B.size1()) {
        throw std::domain_error("The two Matrix have a diffentent number of rows");
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
    find_first_not_null_in_col(const Matrix<R> &A,
                               const size_type &col)
    {
      size_type i=0;

      while ( (i<A.size1()) && (A(i,col)==0) ) {
        ++i;
      }
      return i;

    }


    template <class R>
    Matrix<R>
    remove_null_columns_but_one(const Matrix<R> &A)
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

      Matrix<R> new_A(rows,std::max((size_t)1,directions));

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
    remove_null_columns(const Matrix<R> &A,
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
    Matrix<R>
    compute_space(const Matrix<R> &SA,
                  array<size_type> &row,const array<size_type> &col)
    {
      size_type cols=col.size(), rows=row.size();

      assert(cols>=rows);

      size_type SA_cols=SA.size2(), A_rows=SA_cols-rows;
      size_type i,j,k,j2;

      Matrix<R> A(A_rows,SA_cols);

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
