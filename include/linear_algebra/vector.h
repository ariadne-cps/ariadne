/***************************************************************************
 *            vector.h
 *
 *  Mon May  3 12:31:15 2004
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
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
 
/*! \file vector.h
 *  \brief Vector types and vector operations.
  */

#ifndef _ARIADNE_VECTOR_H
#define _ARIADNE_VECTOR_H 

#include <iosfwd>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../declarations.h"
#include "../numeric/integer.h"



namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A vector over \a R. */
    template<typename R>
    class Vector
      : public boost::numeric::ublas::vector<R> 
    {
      typedef boost::numeric::ublas::vector<R> _Base;
     public:
      Vector() : _Base() { }
      Vector(const size_type& n) : _Base(n) { }
      template<typename E> Vector(const boost::numeric::ublas::vector_expression<E>& v) : _Base(v()) { }
      explicit Vector(const std::string& s); 
    };
    
    
    template <typename R>
    inline Vector<R> zero_vector(dimension_type n) {
      Vector<R> v(n);
      for (dimension_type i=0; i<n; ++i) {
        v(i)=0;
      }
      return v;
    }
    
    template <typename R>
    inline Vector<R> ones_vector(dimension_type n) {
      Vector<R> v(n);
      for (dimension_type i=0; i<n; ++i) {
        v(i)=1;
      }
      return v;
    }
    
    template <typename R>
    inline Vector<R> unit_vector(dimension_type n, dimension_type i) {
      Vector<R> v(n);
      v(i)=1;
      return v;
    }
    
    template <typename R>
    inline R inner_product(const Vector<R> u, const Vector<R>& v) {
      assert(u.size()==v.size());
      R result=0;
      for(size_type i=0; i!=u.size(); ++i) {
        result+=u(i)*v(i);
      }
      return result;
    }

    template <typename R>
    inline bool linear_multiple(const Vector<R> u, const Vector<R>& v) {
      assert(u.size()==v.size());
      R multiple=0;
      for (dimension_type i=0; i<u.dimension(); ++i) {
        if(u(i)==0 && multiple!=0 && v(i)!=0) {
          return false;
        }
        if(multiple==0) {
          multiple=v(i)/u(i);
        }
        else {
          if(multiple*u(i)!=v(i)) {
            return false;
          }
        }
      }
      return true;
    }
    
    template<typename R>
    inline
    R
    norm(const Vector<R>& v) {
      R norm=R(0);
      for (size_type i=0; i<v.size(); ++i) {
        norm=std::max(norm,R(abs(v(i))));
      }
      return norm;
    }
    
    template <typename R>
    inline Integer common_denominator(const Vector<R>& b) 
    {
      Integer denom=1;
      for (dimension_type i=0; i< b.size(); ++i) {
        denom=lcm( denom, denominator(b(i)) );
      }
      return denom;
    }

    template <typename R>
    inline R norm_infinite(const Vector<R>& b) 
    {
      R norm=0.0;
      for (size_t i=0; i< b.size(); i++) {
        if (abs(b[i])>norm) norm=abs(b[i]);
      }
      return norm;
    }
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Vector<R>& v);
    
    template <typename R>
    std::istream&
    operator>>(std::istream& is, Vector<R>& v);

  }
}


#endif /* _ARIADNE_VECTOR_H */
