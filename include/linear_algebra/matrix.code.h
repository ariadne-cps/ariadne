/***************************************************************************
 *            matrix.code.h
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

    template<class R>
    Matrix<R>::Matrix(const std::string& s)
    {
        std::stringstream ss(s);
        ss >> *this;
    }

    template<class R>
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

    template<class R>
    bool
    Matrix<R>::operator!=(const Matrix<R>& A) const
    {
      return !(*this==A);
    }

    template<class R>
    typename Matrix<R>::F
    Matrix<R>::norm() const
    {
      //std::cerr << "Matrix<" << name<R>() << ">::norm() const\n";
      const Matrix<R>& A=*this;
      F result=0;
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        F row_sum=0;
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          row_sum+=abs(A(i,j));
        }
        result=Numeric::max(result,row_sum);
      }
      return result;
    }

    template<class R>
    typename Matrix<R>::F
    Matrix<R>::log_norm() const
    {
      const Matrix<R>& A=*this;
      F result=0;
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        F row_sum=A(i,i);
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          if(i!=j) {
            row_sum+=abs(A(i,j));
          }
        }
        result=max(result,row_sum);
      }
      return result;
    }

    template<class R>
    Matrix<R>
    Matrix<R>::transpose() const
    {
      return Matrix<R>(this->number_of_columns(),this->number_of_rows(),const_cast<R*>(this->begin()),
                       this->column_increment(),this->row_increment());
    }


    template<class R>
    bool
    Matrix<R>::singular() const {
      return LUMatrix<F>(*this).singular();
    }

    template<class R>
    typename Matrix<R>::F
    Matrix<R>::determinant() const
    {
      return LUMatrix<F>(*this).determinant();
    }

    template<class R>
    Matrix<R>
    Matrix<R>::concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2) {
      if(!(A1.number_of_rows()==A2.number_of_rows())) { throw IncompatibleSizes(__PRETTY_FUNCTION__); }
      LinearAlgebra::Matrix<R> result(A1.number_of_rows(),A1.number_of_columns()+A2.number_of_columns());
      for(size_type i=0; i!=result.number_of_rows(); ++i) {
        for(size_type j=0; j!=A1.number_of_columns(); ++j) {
          result(i,j)=A1(i,j);
        }
        for(size_type j=0; j!=A2.number_of_columns(); ++j) {
          result(i,j+A1.number_of_columns())=A2(i,j);
        }
      }
      return result;
    }



    template<class R>
    Matrix<typename Matrix<R>::F>
    Matrix<R>::inverse() const
    {
      return LUMatrix<F>(*this).inverse();
    }


    template<class R>
    Vector<typename Matrix<R>::F>
    Matrix<R>::solve(const Vector<R>& v) const
    {
      return this->inverse()*v;
    }


    template<class R>
    std::ostream&
    Matrix<R>::write(std::ostream& os) const
    {
      const Matrix<R>& A=*this;
      os << "[";
      for(uint i=0; i!=A.number_of_rows(); ++i) {
        for(uint j=0; j!=A.number_of_columns(); ++j) {
          os << (j==0 ? (i==0 ? " " : "; ") : ",");
          //os << Ariadne::convert_to<double>(A(i,j));
          os << A(i,j);
        }
      }
      os << " ]";
      return os;
    }


    template<class R>
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
          for(size_type i=0; i!=A.number_of_rows(); ++i) {
            check_size(v[i],A.number_of_columns(),__PRETTY_FUNCTION__);
            for(size_type j=0; j!=A.number_of_columns(); ++j) {
              A(i,j)=v[i][j];
            }
          }
        }
      }
      else {
        std::cerr << "c=" << c << std::endl;
        throw invalid_input("Matrix<R>");
      }
      return is;
    }





    // FIXME:
    // The following functions are only used by the Zonotope class.
    // Since they are not very general, we should probably remove 
    // many of these operations.
    // Also note that many algorithms can probably be better implemented
    // by using LU or QR factorizations.

    template<class R>
    bool
    independent_rows(Matrix<R> A)
    {
     const size_type rows = A.number_of_rows(),cols = A.number_of_columns();
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


    template<class R>
    bool
    equivalent_columns(const Matrix<R> &A, const size_type &A_col,
                       const Matrix<R> &B, const size_type &B_col)
    {
      if (A.number_of_rows()!=B.number_of_rows()) {
        throw std::domain_error("The two Matrix have a diffentent number of rows");
      }
      for (size_type i=0; i<A.number_of_rows(); i++) {
        if (A(i,A_col)!=B(i,B_col)) {
          return false;
        }
      }
      return true;
    }


    template<class R>
    size_type
    find_first_not_null_in_col(const Matrix<R> &A,
                               const size_type &col)
    {
      size_type i=0;

      while ( (i<A.number_of_rows()) && (A(i,col)==0) ) {
        ++i;
      }
      return i;

    }


    template<class R>
    Matrix<R>
    remove_null_columns_but_one(const Matrix<R> &A)
    {
      size_type directions=0;
      size_type rows=A.number_of_rows();
      size_type cols=A.number_of_columns();
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


    template<class R>
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


    template<class R>
    Matrix<R>
    compute_space(const Matrix<R> &SA,
                  array<size_type> &row,const array<size_type> &col)
    {
      size_type cols=col.size(), rows=row.size();

      if(cols<rows) { throw IncompatibleSizes(__PRETTY_FUNCTION__); }

      size_type SA_cols=SA.number_of_columns(), A_rows=SA_cols-rows;
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
