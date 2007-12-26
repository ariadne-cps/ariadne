/***************************************************************************
 *            vector_slice.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file vector_slice.h
 *  \brief Vector slices
  */

#ifndef ARIADNE_VECTOR_SLICE_H
#define ARIADNE_VECTOR_SLICE_H 

#include <iosfwd>
#include <algorithm>

#include "base/types.h"
#include "base/array.h"
#include "numeric/declarations.h"
#include "numeric/integer.h"
#include "numeric/arithmetic.h"
#include "numeric/interval.h"

#include "linear_algebra/exceptions.h"
#include "linear_algebra/vector_expression.h"

namespace Ariadne {
  namespace LinearAlgebra {
    
    class Slice;
    template<class R> class Vector;
    template<class R> class VectorSlice;
    template<class R> class VectorSlice<const R>;


    /*! \ingroup LinearAlgebra
     *  \brief A slice through an array or vector, with equally spaces increments. 
     */
    template<class R>
    class VectorSlice
      : public VectorExpression< VectorSlice<R> >
    {
     public:
      typedef R value_type;
     
      explicit VectorSlice(const size_type& size, R* begin, const size_type& increment=1u)
        : _size(size), _begin(begin), _increment(increment) { }
      VectorSlice(Vector<R>& v)
        : _size(v.size()), _begin(v.data().begin()), _increment(1u) { }
      size_type size() const { return this->_size; }
      size_type increment() const { return this->_increment; };
      R* begin() { return this->_begin; }
      const R* begin() const { return this->_begin; } 
      R& operator() (const size_type& i) { return this->_begin[i*this->_increment]; }
      R& operator[] (const size_type& i) { return this->_begin[i*this->_increment]; }
      const R& operator() (const size_type& i) const { return this->_begin[i*this->_increment]; };
      const R& operator[] (const size_type& i) const { return this->_begin[i*this->_increment]; }
     
      template<class E> VectorSlice<R>& operator=(const VectorExpression< E >& v);
     private:
      size_type _size;
      R* _begin;
      size_type _increment;
    };
    

    template<class R>
    class VectorSlice<const R>
      : public VectorExpression< VectorSlice<const R> >
    {
     public:
      typedef R value_type;
      explicit VectorSlice(const size_type& size, const R* begin, const size_type& increment=1u)
        : _size(size), _begin(begin), _increment(increment) { }
      VectorSlice(const Vector<R>& v)
        : _size(v.size()), _begin(v.data().begin()), _increment(1u) { }
      size_type size() const { return this->_size; }
      size_type increment() const { return this->_increment; };
      const R* begin() const { return this->_begin; } 
      const R& operator() (const size_type& i) const { return this->_begin[i*this->_increment]; };
      const R& operator[] (const size_type& i) const { return this->_begin[i*this->_increment]; }
     private:
      // Disable assignment operators
      VectorSlice<const R>& operator=(const Vector<R>&);
      VectorSlice<const R>& operator=(const VectorSlice<const R>&);
     private:
      size_type _size;
      const R* _begin;
      size_type _increment;
    };
    
    template<class R> template<class E> inline
    VectorSlice<R>&
    VectorSlice<R>::operator=(const VectorExpression< E >& v) 
    {
      const E& e=v(); 
      ARIADNE_CHECK_EQUAL_SIZES(*this,e,"VectorSlice& VectorSlice::operator=(VectorExpression)");
      for(size_type i=0; i!=e.size(); ++i) {
        this->_begin[i*this->_increment]=e(i); 
      }
      return *this;
    }

    template<class R> std::ostream& operator<<(std::ostream& os, const VectorSlice<R>& v) {
      return os << Vector<R>(v); }
    template<class R> std::ostream& operator<<(std::ostream& os, const VectorSlice<const R>& v) {
      return os << Vector<R>(v); }

  }
}

#endif /* ARIADNE_VECTOR_SLICE_H */
