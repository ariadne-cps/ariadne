/***************************************************************************
 *            interval_matrix.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file interval_matrix.h
 *  \brief Matrices of intervals.
  */

#ifndef _ARIADNE_INTERVAL_MATRIX_H
#define _ARIADNE_INTERVAL_MATRIX_H 

#include <iosfwd>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "../linear_algebra/matrix.h"
#include "../declarations.h"

#include "../numeric/interval.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A matrix of intervals. */
    template<typename R>
    class IntervalMatrix : public Matrix< Interval<R> >
    {
     private:
      typedef Matrix< Interval<R> > _Base;
     public:
      /*! \brief Construct a 0 by 0 interval matrix. */
      explicit IntervalMatrix() : _Base() { }
      /*! \brief Construct an \a r by \a c interval matrix, all of whose entries are zero. */
      explicit IntervalMatrix(const size_type& r, const size_type& c) : _Base(r,c) { 
        Interval<R> z(0); for(size_type i=0; i!=r; ++i) { for(size_type j=0; j!=c; ++j) { (*this)(i,j)=z; } } 
      }
      /*! \brief Construct an \a r by \a c interval matrix, from the array starting at \a ptr. */
      explicit IntervalMatrix(const size_type& r, const size_type& c, const Interval<R>* ptr, 
                              const size_type& ir, const size_type& ic) : _Base(r,c,ptr,ir,ic) { }
      /*! \brief Construct an interval vector centred at \a A, each of whose elements has radius \a r. */
      explicit IntervalMatrix(const Matrix<R>& A, const R& e) : _Base(A) { 
        size_type m=this->number_of_rows(); size_type n=this->number_of_columns();
        R r=e/n; Interval<R> ir(-r,r); 
        for(size_type i=0; i!=m; ++i) { for(size_type j=0; j!=n; ++j) { (*this)(i,j)+=ir; } } 
      }
      
      /*! \brief Convert from a matrix. */
      IntervalMatrix(const Matrix<R>& A) : _Base(A) { }
      
      /* Convert from a matrix expression. */
      template<typename E> IntervalMatrix(const boost::numeric::ublas::matrix_expression<E>& A) : _Base(A()) { }

      /*! \brief Construct from a string literal of the form 
       *  "[[l11,u11],[l12,u12],...,[l1n,u1n]; [l21,u21],[l22,u22],...,[l2n,u2n];...;[lm1,um1],[lm2,um2],...,[lmn,umn]]". */
      explicit IntervalMatrix(const std::string& s) : _Base(s) { }
      
      /*! \brief An \a r by \a c interval matrix, all of whose entries are zero. */
      static IntervalMatrix<R> zero(const size_type r, const size_type c);
      /*! \brief The \a n by \a n identity interval matrix. */
      static IntervalMatrix<R> identity(const size_type n);
      
      /*! \brief The centre of the interval matrix. */
      Matrix<R> centre() const;
      /*! \brief The sums of the radii in each row. */
      IntervalVector<R> radius_row_sum() const;
      /*! \brief The maximum of the sums of the radii in each row. */
      R radius_norm() const;
      
      /*! \brief The maximum possible norm. */
      R upper_norm() const;
      /*! \brief The maximum possible logarithmic norm. */
      R upper_log_norm() const;
    };



    /*! \brief A matrix \f$A\f$ such that for all zonotopes \f$Z\f$, \f$AZ\subset \overline{\underline{A}}\f$. */
    template<typename R>
    Matrix<R>
    over_approximation(const IntervalMatrix<R>& A); 
        
    /*! \brief An interval matrix exponential. */
    template<typename R>
    IntervalMatrix<R>
    exp(const IntervalMatrix<R>& A); 
        
    /*! \brief The interval matrix inverse. */
    template<typename R>
    IntervalMatrix<R>
    approximate(const Matrix<typename numerical_traits<R>::field_extension_type>& A,const R& e); 
        
    
    template <typename R>
    IntervalVector<R> 
    prod(const Matrix<R>& A, const IntervalVector<R>& v);
        
    template <typename R>
    IntervalVector<R> 
    prod(const IntervalMatrix<R>& A, const Vector<R>& v);
        
    template <typename R>
    IntervalVector<R> 
    prod(const IntervalMatrix<R>& A, const IntervalVector<R>& v);
        
    template <typename R>
    IntervalMatrix<R> 
    prod(const IntervalMatrix<R>& A, const Matrix<R>& B);
      
    template <typename R>
    IntervalMatrix<R> 
    prod(const Matrix<R>& A, const IntervalMatrix<R>& B);
      
    template <typename R>
    IntervalMatrix<R> 
    prod(const IntervalMatrix<R>& A, const IntervalMatrix<R>& B);
      
    template<typename R>
    IntervalMatrix<R>
    fprod(const Matrix<typename numerical_traits<R>::field_extension_type>& A, 
         const IntervalMatrix<R>& B);
    

    template<typename R>
    inline
    IntervalVector<R>
    operator*(const Matrix<R>& A, const IntervalVector<R>& B) {
      return prod(A,B);
    }
      
    template<typename R>
    inline
    IntervalVector<R>
    operator*(const IntervalMatrix<R>& A, const Vector<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    IntervalVector<R>
    operator*(const IntervalMatrix<R>& A, const IntervalVector<R>& B) {
      return prod(A,B);
    }
    
    
    template<typename R>
    inline
    IntervalMatrix<R>
    operator*(const IntervalMatrix<R>& A, const Matrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    IntervalMatrix<R>
    operator*(const Matrix<R>& A, const IntervalMatrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    IntervalMatrix<R>
    operator*(const IntervalMatrix<R>& A, const IntervalMatrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    IntervalMatrix<R>
    operator*(const Matrix<typename numerical_traits<R>::field_extension_type>& A, 
              const IntervalMatrix<R>& B) {
      return fprod(A,B);
    }

        
    /*! \brief The range of supremum norms of matrices in the interval Matrix. */
    template<typename R>
    inline
    Interval<R>
    norm(const IntervalMatrix<R>& A)
    { 
      return A.norm();
    }
        
    /*! \brief The maximum norm in the interval matrix.
     */
    template<typename R>
    inline
    R
    upper_norm(const IntervalMatrix<R>& A)
    {    
      return A.upper_norm();
    }
    
    /*! \brief The maximum logarithmic norm in the interval matrix.
     * 
     *  The logarithmic norm is defined as \f$\max_{i} A_{ii}+\sum_{j!=i}|A_{ij}|\f$.
     */
    template<typename R>
    inline
    R
    upper_log_norm(const IntervalMatrix<R>& A)
    {
      return A.upper_log_norm();
    }
        
    /*! \brief The interval matrix inverse. */
    template<typename R>
    inline
    IntervalMatrix<R>
    inverse(const IntervalMatrix<R>& A)
    {
      return A.inverse(); 
    }
    
    template <typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const IntervalMatrix<R>& A)
    {
      return A.write(os);
    }
    
    template <typename R>
    inline
    std::istream&
    operator>>(std::istream& is, IntervalMatrix<R>& A)
    {
      return A.read(is);
    }
    
    
  }
}

#endif /* _ARIADNE_INTERVAL_MATRIX_H */
