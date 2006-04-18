/***************************************************************************
 *            tensor.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file tensor.h
 *  \brief Tensors.
 */

#ifndef _ARIADNE_TENSOR_H
#define _ARIADNE_TENSOR_H

#include <iosfwd>

#include "../declarations.h"
#include "../base/array.h"


namespace Ariadne {
  namespace LinearAlgebra {
    
    
    /*! \brief A Tensor representing a second derivative. */
    template<typename R> 
    class Tensor {
     public:
      Tensor(size_type m, size_type n1, size_type n2)
        : _sizes(3), _elements(m*n1*n2,R(0)) 
      { _sizes[0]=m; _sizes[1]=n1; _sizes[2]=n2; }
      
      Tensor(const Tensor<R>& other) : _sizes(other._sizes), _elements(other._elements) { }

      Tensor<R>& operator=(const Tensor<R>& other) {
        if(this!=&other) {
          this->_sizes=other._sizes;
          this->_elements=other._elements;
        }
        return *this;
      }        
      
      size_type size(const size_type& i) const { return _sizes[i]; }
      
      R operator() (const size_type& i, const size_type& j1, const size_type& j2) const {
        return this->_elements[(i*_sizes[1]+j1)*_sizes[2]+j2];
      }
      
      R& operator() (const size_type& i, const size_type& j1, const size_type& j2) {
        return this->_elements[(i*_sizes[1]+j1)*_sizes[2]+j2];
      }
     private:
      array<size_type> _sizes;
      array<R> _elements;
    };
  
    template<typename R>
    Vector<R>
    product(const Tensor<R>& T, const Vector<R>& v1, const Vector<R>& v2);
    
    template<typename R>
    Matrix<R>
    product(const Tensor<R>& T, const Vector<R>& v1);
    
    template<typename R>
    Tensor<R>
    product(const Tensor<R>& T, const Matrix<R>& A1);

    template<typename R>
    std::ostream&
    operator<< (std::ostream& os, const Tensor<R>& T); 

    template<typename R>
    inline 
    Matrix<R> 
    operator*(const Tensor<R>& T, const Vector<R>& v)
    { 
      return product(T,v); 
    }
    
  }
}

#endif /* _ARIADNE_TENSOR_H */
