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
    class interval_matrix : public boost::numeric::ublas::matrix< Interval<R> >
    {
     private:
      typedef boost::numeric::ublas::matrix< Interval<R> > Base;
     public:
      interval_matrix() : Base() { }
      interval_matrix(const size_type& r, const size_type& c) : Base(r,c) { }
      template<typename E> interval_matrix(const boost::numeric::ublas::matrix_expression<E>& A) : Base(A()) { }
      interval_matrix(const matrix<R>& A, const R& r);
      
      matrix<R> centre() const;
      R radius() const;
      
      Interval<R> norm() const;
      R upper_norm() const;
      R upper_log_norm() const;
      
      interval_matrix<R> inverse() const;
    };

    template <typename R>
    std::ostream& 
    operator<<(std::ostream& os, const interval_matrix<R>& A);
        
    /*! \brief A matrix \f$A\f$ such that for all zonotopes \f$Z\f$, \f$AZ\subset \overline{\underline{A}}\f$. */
    template<typename R>
    matrix<R>
    over_approximation(const interval_matrix<R>& A); 
        
    /*! \brief An interval matrix exponential. */
    template<typename R>
    interval_matrix<R>
    exp(const interval_matrix<R>& A); 
        
    /*! \brief The interval matrix inverse. */
    template<typename R>
    interval_matrix<R>
    approximate(const matrix<typename numerical_traits<R>::field_extension_type>& A,const R& e); 
        
    
    template <typename R>
    interval_vector<R> 
    prod(const matrix<R>& A, const interval_vector<R>& v);
        
    template <typename R>
    interval_vector<R> 
    prod(const interval_matrix<R>& A, const vector<R>& v);
        
    template <typename R>
    interval_vector<R> 
    prod(const interval_matrix<R>& A, const interval_vector<R>& v);
        
    template <typename R>
    interval_matrix<R> 
    prod(const interval_matrix<R>& A, const matrix<R>& B);
      
    template <typename R>
    interval_matrix<R> 
    prod(const matrix<R>& A, const interval_matrix<R>& B);
      
    template <typename R>
    interval_matrix<R> 
    prod(const interval_matrix<R>& A, const interval_matrix<R>& B);
      
    template<typename R>
    interval_matrix<R>
    fprod(const matrix<typename numerical_traits<R>::field_extension_type>& A, 
         const interval_matrix<R>& B);
    

    template<typename R>
    inline
    interval_vector<R>
    operator*(const matrix<R>& A, const interval_vector<R>& B) {
      return prod(A,B);
    }
      
    template<typename R>
    inline
    interval_vector<R>
    operator*(const interval_matrix<R>& A, const vector<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    interval_vector<R>
    operator*(const interval_matrix<R>& A, const interval_vector<R>& B) {
      return prod(A,B);
    }
    
    
    template<typename R>
    inline
    interval_matrix<R>
    operator*(const interval_matrix<R>& A, const matrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    interval_matrix<R>
    operator*(const matrix<R>& A, const interval_matrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    interval_matrix<R>
    operator*(const interval_matrix<R>& A, const interval_matrix<R>& B) {
      return prod(A,B);
    }
    
    template<typename R>
    inline
    interval_matrix<R>
    operator*(const matrix<typename numerical_traits<R>::field_extension_type>& A, 
              const interval_matrix<R>& B) {
      return fprod(A,B);
    }

    /*! \brief The range of supremum norms of matrices in the interval matrix. */
    template<typename R>
    inline
    Interval<R>
    norm(const interval_matrix<R>& A)
    { 
      return A.norm();
    }
        
    /*! \brief The maximum norm in the interval matrix.
     */
    template<typename R>
    inline
    R
    upper_norm(const interval_matrix<R>& A)
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
    upper_log_norm(const interval_matrix<R>& A)
    {
      return A.upper_log_norm();
    }
        
    /*! \brief The interval matrix inverse. */
    template<typename R>
    inline
    interval_matrix<R>
    inverse(const interval_matrix<R>& A)
    {
      return A.inverse(); 
    }
    
  }
}

#endif /* _ARIADNE_INTERVAL_MATRIX_H */
