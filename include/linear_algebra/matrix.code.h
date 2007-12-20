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

#include <algorithm>

#include "base/array.h"

#include "numeric/declarations.h"
#include "numeric/integer.h"
#include "numeric/arithmetic.h"
#include "numeric/approximate_float.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/lu_matrix.h"
#include "linear_algebra/qr_matrix.h"

namespace Ariadne {


namespace {

template<class R>
void
instantiate_matrix_approx() 
{
  typedef typename Numeric::traits<R>::number_type X;
  LinearAlgebra::Matrix<X>* nA=0;
  LinearAlgebra::Vector<X>* nv=0;

  LinearAlgebra::mul_approx(*nA,*nv);
  LinearAlgebra::mul_approx(*nA,*nA);
  LinearAlgebra::solve_approx(*nA,*nv);
  LinearAlgebra::inverse_approx(*nA);
  LinearAlgebra::qr_approx(*nA);
}

template<>
void
instantiate_matrix_approx<Numeric::Rational>() 
{
}

template<>
void
instantiate_matrix_approx< Numeric::Interval<Numeric::Rational> >() 
{
}

} // namespace


template<class R>
void
LinearAlgebra::Matrix<R>::instantiate()
{
  typedef typename Numeric::traits<R>::number_type X;
  typedef typename Numeric::traits<R>::interval_type I;

  Matrix<R>* A=0;
  Vector<R>* v=0;
  Vector<X>* nv=0;
  Vector<I>* iv=0;

  sup_norm(*A);
  log_norm(*A);
  singular(*A);
  determinant(*A);
  inverse(*A);
  solve(*A,*nv);
  solve(*A,*iv);
  transpose(*A);
  direct_sum(*A,*A);
  concatenate(*A,*A);
  concatenate_rows(*A,*A);
  concatenate_rows(*A,*v);
  concatenate_columns(*A,*A);
  concatenate_columns(*A,*v);

  instantiate_matrix_approx<R>();
}

template<class R>
LinearAlgebra::Matrix<R>::Matrix(const std::string& s)
{
  std::stringstream ss(s);
  ss >> *this;
}
    

template<class R>
typename Numeric::traits<R>::arithmetic_type
LinearAlgebra::sup_norm(const Matrix<R>& A)
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
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
LinearAlgebra::log_norm(const Matrix<R>& A) 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
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
LinearAlgebra::Matrix<R>
LinearAlgebra::transpose(const Matrix<R>& A)
{
  return Matrix<R>(A.number_of_columns(),A.number_of_rows(),const_cast<R*>(A.begin()),
                   A.column_increment(),A.row_increment());
}


template<class R>
bool
LinearAlgebra::singular(const Matrix<R>& A) 
{
  return LUMatrix<typename Numeric::traits<R>::arithmetic_type>(A).singular();
}


template<class R>
typename Numeric::traits<R>::arithmetic_type
LinearAlgebra::determinant(const Matrix<R>& A) 
{
  return LUMatrix<typename Numeric::traits<R>::arithmetic_type>(A).determinant();
}


template<class R>
LinearAlgebra::Matrix<R>
LinearAlgebra::direct_sum(const Matrix<R>& A1, const Matrix<R>& A2) {
  return concatenate(A1,A2);
}


template<class R>
LinearAlgebra::Matrix<R>
LinearAlgebra::concatenate(const Matrix<R>& A1, const Matrix<R>& A2) {
  LinearAlgebra::Matrix<R> result(A1.number_of_rows()+A2.number_of_rows(),A1.number_of_columns()+A2.number_of_columns());
  for(size_type j=0; j!=A1.number_of_columns(); ++j) {
    for(size_type i=0; i!=A1.number_of_rows(); ++i) {
      result(i,j)=A1(i,j);
    }
  }
  for(size_type j=0; j!=A2.number_of_columns(); ++j) {
    for(size_type i=0; i!=A2.number_of_rows(); ++i) {
      result(i+A1.number_of_rows(),j+A1.number_of_columns())=A2(i,j);
    }
  }
  return result;
}


template<class R>
LinearAlgebra::Matrix<R>
LinearAlgebra::concatenate_rows(const Matrix<R>& A, const Vector<R>& v) {
  if(!(A.number_of_columns()==v.size())) { 
    ARIADNE_THROW(IncompatibleSizes,"Matrix concatenate_rows(Matrix,Matrix)","A="<<A<<", v="<<v); 
  }
  LinearAlgebra::Matrix<R> result(A.number_of_rows()+1u,A.number_of_columns());
  for(size_type j=0; j!=A.number_of_columns(); ++j) {
    for(size_type i=0; i!=A.number_of_rows(); ++i) {
      result(i,j)=A(i,j);
    }
    result(A.number_of_rows(),j)=v(j);
  }
  return result;
}


template<class R>
LinearAlgebra::Matrix<R>
LinearAlgebra::concatenate_rows(const Matrix<R>& A1, const Matrix<R>& A2) {
  if(!(A1.number_of_columns()==A2.number_of_columns())) { 
    ARIADNE_THROW(IncompatibleSizes,"Matrix concatenate_rows(Matrix,Matrix)","A1="<<A1<<", A2="<<A2); 
  }
  LinearAlgebra::Matrix<R> result(A1.number_of_rows()+A2.number_of_rows(),A1.number_of_columns());
  for(size_type j=0; j!=result.number_of_columns(); ++j) {
    for(size_type i=0; i!=A1.number_of_rows(); ++i) {
      result(i,j)=A1(i,j);
    }
    for(size_type i=0; i!=A2.number_of_rows(); ++i) {
      result(i+A1.number_of_rows(),j)=A2(i,j);
    }
  }
  return result;
}


template<class R>
LinearAlgebra::Matrix<R>
LinearAlgebra::concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2) {
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
LinearAlgebra::Matrix<R>
LinearAlgebra::concatenate_columns(const Matrix<R>& A, const Vector<R>& v) {
  if(!(A.number_of_rows()==v.size())) {
    ARIADNE_THROW(IncompatibleSizes,"Matrix concatenate_columns(Matrix,Vector)","A="<<A<<", v="<<v); 
  }
  LinearAlgebra::Matrix<R> result(A.number_of_rows(),A.number_of_columns()+1u);
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(i,j)=A(i,j);
    }
    result(i,A.number_of_columns())=v(i);
  }
  return result;
}



template<class R>
LinearAlgebra::Matrix<typename Numeric::traits<R>::arithmetic_type>
LinearAlgebra::inverse(const Matrix<R>& A) 
{
  return LUMatrix<typename Numeric::traits<R>::arithmetic_type>(A).inverse();
}



template<class T>
LinearAlgebra::Vector< Numeric::Float<T> >
LinearAlgebra::mul_approx(const Matrix< Numeric::Float<T> >& A,
                          const Vector< Numeric::Float<T> >& b)
{
  typedef Numeric::Float<T> F;
  typedef Numeric::ApproximateFloat<T> AF;
  Vector<F> r;
  reinterpret_cast<Vector<AF>&>(r)=reinterpret_cast<Matrix<AF>const&>(A) * reinterpret_cast<Vector<AF>const&>(b);
  return r;
}

template<class T>
LinearAlgebra::Matrix< Numeric::Float<T> >
LinearAlgebra::mul_approx(const Matrix< Numeric::Float<T> >& A,
                          const Matrix< Numeric::Float<T> >& B)
{
  typedef Numeric::Float<T> F;
  typedef Numeric::ApproximateFloat<T> AF;
  Matrix<F> C;
  reinterpret_cast<Vector<AF>&>(C)=reinterpret_cast<Matrix<AF>const&>(A) * reinterpret_cast<Vector<AF>const&>(A);
  return C;
}



template<class T>
LinearAlgebra::Vector< Numeric::Float<T> >
LinearAlgebra::solve_approx(const Matrix< Numeric::Float<T> >& A,
                            const Vector< Numeric::Float<T> >& b)
{
  typedef Numeric::Float<T> F;
  typedef Numeric::ApproximateFloat<T> AF;
  Vector<F> r;
  reinterpret_cast<Vector<AF>&>(r)=solve(reinterpret_cast<Matrix<AF>const&>(A),reinterpret_cast<Vector<AF>const&>(b));
  return r;
}

template<class R1, class R2>
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::solve(const Matrix<R1>& A, const Vector<R2>& v) 
{
  return inverse(A)*v;
}


template<class T>
LinearAlgebra::Matrix< Numeric::Float<T> >
LinearAlgebra::inverse_approx(const Matrix< Numeric::Float<T> >& A)
{
  typedef Numeric::Float<T> F;
  typedef Numeric::ApproximateFloat<T> AF;
  Matrix<F> r;
  reinterpret_cast<Matrix<AF>&>(r)=inverse(reinterpret_cast<Matrix<AF>const&>(A));
  return r;
}

template<class T>
std::pair< LinearAlgebra::Matrix< Numeric::Float<T> >, 
           LinearAlgebra::Matrix< Numeric::Float<T> > >
LinearAlgebra::qr_approx(const Matrix< Numeric::Float<T> >& A)
{
  typedef Numeric::Float<T> F;
  typedef Numeric::ApproximateFloat<T> AF;
  QRMatrix<AF> qr(reinterpret_cast<Matrix<AF>const&>(A));
  Matrix<F> Q;
  Matrix<F> R;
  Matrix<AF>& AQ=reinterpret_cast<Matrix<AF>&>(Q);
  Matrix<AF>& AR=reinterpret_cast<Matrix<AF>&>(R);
  AQ=qr.Q();
  AR=qr.R();
  return std::pair< Matrix<F>, Matrix<F> >(Q,R);
}

template<class R>
LinearAlgebra::Matrix<typename Numeric::traits<R>::arithmetic_type>
LinearAlgebra::schulz_inverse(const Matrix<R>& A) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
std::ostream&
LinearAlgebra::Matrix<R>::write(std::ostream& os) const
{
  const Matrix<R>& A=*this;
  os<<std::fixed<<std::setprecision(6);
  os << "[";
  for(uint i=0; i!=A.number_of_rows(); ++i) {
    for(uint j=0; j!=A.number_of_columns(); ++j) {
      os << (j==0 ? (i==0 ? "" : "; ") : ",");
      //os << Ariadne::convert_to<double>(A(i,j));
      os << A(i,j);
    }
  }
  os << " ]\n";
  return os;
}


template<class R>
std::istream&
LinearAlgebra::Matrix<R>::read(std::istream& is)
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
