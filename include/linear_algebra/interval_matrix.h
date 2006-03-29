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

#include "../base/basic_type.h"
#include "../base/numerical_type.h"
#include "../base/interval.h"

namespace boost {
  namespace numeric {
    namespace ublas {
      
      using Ariadne::Interval;
      using Ariadne::size_type;
      
      template <typename Real>
      vector< Interval<Real> > 
      iprod(const matrix<Real>& A, const vector< Interval<Real> >& B) {
        vector< Interval<Real> > result(A.size1());
        for (size_type i=0; i!=result.size(); ++i) {
          result(i)=Interval<Real>(0);
          for (size_type j=0; j!=B.size(); ++j) {
            result(i)+=A(i,j)*B(j);
          }
        }
        return result;
      }
        
      template <typename Real>
      matrix< Interval<Real> > 
      iprod(const matrix< Interval<Real> >& A, const matrix<Real>& B) {
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
      iprod(const matrix<Real>& A, const matrix< Interval<Real> >& B) {
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
      
      template<typename R>
      inline
      vector< Interval<R> >
      operator*(const matrix<R>& A, const vector< Interval<R> >& B) {
        return iprod(A,B);
      }
      
      
      template<typename R>
      inline
      matrix< Interval<R> >
      operator*(const matrix< Interval<R> >& A, const matrix<R>& B) {
        return iprod(A,B);
      }
      
      template<typename R>
      inline
      matrix< Interval<R> >
      operator*(const matrix<R>& A, const matrix< Interval<R> >& B) {
        return iprod(A,B);
      }
   

    
    }
  }
}



#endif /* _ARIADNE_INTERVAL_MATRIX_H */
