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
#include "../numeric/arithmetic.h"



namespace Ariadne {
  namespace LinearAlgebra {
    
    /*! \brief A vector over \a R. */
    template<typename R>
    class Vector
      : public boost::numeric::ublas::vector<R> 
    {
      typedef boost::numeric::ublas::vector<R> _Base;
     public:
      /*! \brief Construct a vector of size 0. */
      explicit Vector() : _Base() { }
      /*! \brief Construct the zero vector of size \a n. */
      explicit Vector(const size_type& n) : _Base(n) { }
      /*! \brief Construct a vector of size \a n from the array beginning at \a ptr. */
      explicit Vector(const size_type& n, const R* ptr) : _Base(n) { 
        for(size_type i=0; i!=n; ++i) { (*this)(i)=ptr[i]; } }

      /* Convert from a vector expression. */
      template<typename E> Vector(const boost::numeric::ublas::vector_expression<E>& v) : _Base(v) { }
        
      /*! \brief Construct from a string literal of the form 
       *  "[v1,v2,...,vn]", where the vi are numerical literals for the vector elements. 
       */
      explicit Vector(const std::string& s); 

#ifdef DOXYGEN
      /*! \brief Copy constructor. */
      Vector(const Vector<R>& v);
      /*! \brief Copy assignment operator. */
      Vector<R>& operator=(const Vector<R>& v);
#endif
      
      /*! \brief The equality operator. */
      bool operator==(const Vector<R>& v) const {
        if(this->size()!=v.size()) { return false; }
        for(size_type i=0; i!=this->size(); ++i) { if((*this)(i)!=v(i)) { return false; } }
        return true; }
        
      /*! \brief The inequality operator. */
      bool operator!=(const Vector<R>& v) const {
        return !(*this==v); }
      
      /*! \brief Returns the zero vector of size \a n. */
      static Vector<R> zero(const size_type& n) {
        Vector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(0); } return v; }
      /*! \brief Returns the vector of size \a n whose elements are all 1. */
      static Vector<R> one(const size_type& n) {
        Vector<R> v(n); for (dimension_type i=0; i<n; ++i) { v(i)=R(1); } return v; }
      /*! \brief Returns the unit vector of size \a n whose \a i th element is 1. */
      static Vector<R> unit(dimension_type n, dimension_type i) {
        Vector<R> v(n); for (dimension_type k=0; k<n; ++k) { v(k)=R(0); } v(i)=R(1); return v; }
        
#ifdef DOXYGEN
      /*! \brief The number of elements of the vector. */
      size_type size() const;
      /*! \brief A reference to the \a i th element. */
      R& operator() (const size_type& i);
      /*! \brief A constant reference to the \a i th element. */
      const R& operator() (const size_type& i) const;
#endif
        
      /*! \brief The supremum norm. */
      R norm() const {
        R result=0; for (size_type i=0; i<this->size(); ++i) {
          R absi=abs((*this)(i)); result=Numeric::max(result,absi); } return result; }
        
#ifdef DOXYGEN
      /*! \brief The additive inverse of the vector \a v. */
      friend Vector<R> operator-<>(const Vector<R>& v);
      /*! \brief The vector sum of \a v1 and \a v2. */
      friend Vector<R> operator+<>(const Vector<R>& v1, const Vector<R>& v2);
      /*! \brief The vector difference of \a v1 and \a v2. */
      friend Vector<R> operator-<>(const Vector<R>& v1, const Vector<R>& v2);
      /*! \brief The scalar product of \a v by \a s. */
      friend Vector<R> operator*<>(const R& s, const Vector<R>& v);
      /*! \brief The scalar product of \a v by \a s. */
      friend Vector<R> operator*<>(const Vector<R>& v, const R& s);
      /*! \brief The scalar product of \a v by the reciprocal of \a s. */
      friend Vector<R> operator/<>(const Vector<R>& v, const R& s);

      /*! \brief True if \a v1 and \a v2 are parallel. */
      friend bool linear_multiple<>(const Vector<R>& v1, const Vector<R>& v2);
      /*! \brief The inner product of \a v1 and \a v2. */
      friend R inner_product<>(const Vector<R>& v1, const Vector<R>& v2);
      /*! \brief The supremum norm of \a v. */
      friend R LinearAlgebra::norm<>(const Vector<R>& v);
#endif 

      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);
    };
    
    template <typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const Vector<R>& v)
    {
       return v.write(os);
    }
    
    template <typename R>
    inline
    std::istream&
    operator>>(std::istream& is, Vector<R>& v)
    {
       return v.read(is);
    }
    
    template <typename R>
    inline R inner_product(const Vector<R>& u, const Vector<R>& v) 
    {
      assert(u.size()==v.size());
      R result=0;
      for(size_type i=0; i!=u.size(); ++i) {
        result+=u(i)*v(i);
      }
      return result;
    }

    template <typename R>
    inline bool linear_multiple(const Vector<R>& u, const Vector<R>& v) 
    {
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
    inline R norm(const Vector<R>& v) 
    {
      return v.norm(); 
    }
    
  }
}


#endif /* _ARIADNE_VECTOR_H */
