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
    typename Numeric::traits<R>::arithmetic_type
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
    typename Numeric::traits<R>::arithmetic_type
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
      if(!(A1.number_of_rows()==A2.number_of_rows())) { 
        ARIADNE_THROW(IncompatibleSizes,"Matrix concatenate_columns(Matrix,Matrix)","A1="<<A1<<", A2="<<A2); 
      }
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
      os << "[" <<std::endl;
      for(uint i=0; i!=A.number_of_rows(); ++i) {
        for(uint j=0; j!=A.number_of_columns(); ++j) {
          os << (j==0 ? (i==0 ? " " : ";\n ") : ",");
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
            if(v[i].size()!=A.number_of_columns()) {
              ARIADNE_THROW(InvalidInput,"Matrix::read(istream&)","row[0].size()="<<v[0].size()<<", row["<<i<<"].size()="<<v[i].size());
            }
            for(size_type j=0; j!=A.number_of_columns(); ++j) {
              A(i,j)=v[i][j];
            }
          }
        }
      }
      else {
        ARIADNE_THROW(InvalidInput,"Matrix::read(istream&)"," separator c="<<c);
      }
      return is;
    }


  }
}
