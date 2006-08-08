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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "../linear_algebra/matrix.h"
#include "../declarations.h"

#include "../numeric/interval.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A matrix of intervals. */
    template<typename R>
    class IntervalMatrix : public boost::numeric::ublas::matrix< Interval<R> >
    {
     private:
      typedef boost::numeric::ublas::matrix< Interval<R> > Base;
     public:
      IntervalMatrix() : Base() { }
      IntervalMatrix(const size_type& r, const size_type& c) : Base(r,c) { }
      template<typename E> IntervalMatrix(const boost::numeric::ublas::matrix_expression<E>& A) : Base(A()) { }
      IntervalMatrix(const Matrix<R>& A, const R& r);
      
      Matrix<R> centre() const;
      R radius() const;
      
      Interval<R> norm() const;
      R upper_norm() const;
      R upper_log_norm() const;
      /*! \brief Sums of the radii in each row. */
      IntervalVector<R> radius_row_sum() const;
    
      IntervalMatrix<R> inverse() const;
    };

    template <typename R>
    std::ostream& 
    operator<<(std::ostream& os, const IntervalMatrix<R>& A);
        
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
    
  }
}

#endif /* _ARIADNE_INTERVAL_MATRIX_H */
