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
    
    
    /*! \ingroup LinearAlgebra
     *  \brief A tensor representing a second derivative. 
     */
    template<typename R> 
    class Tensor {
     public:
      /*! \brief Construct a rank-three tensor with dimensions \a m, \a n1 and \a n2. */
      Tensor(size_type m, size_type n1, size_type n2)
        : _sizes(3), _strides(3), _elements(m*n1*n2,R(0)) 
      { _sizes[0]=m; _sizes[1]=n1; _sizes[2]=n2; 
        _strides=compute_strides(_sizes); }
      
      /*! \brief Construct a tensor with sizes \a sz from the elements beginning at \a ptr. */
      Tensor(const array<size_type>& sz, const R* ptr) 
       : _sizes(sz), _strides(compute_strides(sz)), _elements(ptr,ptr+this->number_of_elements()) 
      { }
      
      /*! \brief Construct from a string. */
      Tensor(const std::string& str);
      
      /*! \brief Copy constructor. */
      Tensor(const Tensor<R>& t) : _sizes(t._sizes), _elements(t._elements) { }

      /*! \brief Copy assignment operator. */
      Tensor<R>& operator=(const Tensor<R>& t) {
        if(this!=&t) {
          this->_sizes=t._sizes;
          this->_elements=t._elements;
        }
        return *this;
      }        
      
      /*! \brief The rank of the tensor. */
      size_type rank() const { return _sizes.size(); }
      
      /*! \brief The rank of the tensor. */
      size_type number_of_elements() const { 
        return _strides[0]*_sizes[0]; }
        
      /*! \brief The array of sizes in each dimension. */
      const array<size_type>& sizes() const { return _sizes; }
      /*! \brief The size of the \a i th dimension. */
      size_type size(const size_type& i) const { return _sizes[i]; }
      
      /*! \brief A constant reference to the (\a i,\a j1,\a j2)-th element. */
      const R& operator() (const size_type& i, const size_type& j1, const size_type& j2) const {
        return this->_elements[(i*_sizes[1]+j1)*_sizes[2]+j2];
      }
      
      /*! \brief A reference to the (\a i,\a j1,\a j2)-th element. */
      R& operator() (const size_type& i, const size_type& j1, const size_type& j2) {
        return this->_elements[(i*_sizes[1]+j1)*_sizes[2]+j2];
      }
      
      /*! \brief A constant reference to the element indexed by \a ind. */
      const R& operator() (const array<size_type> ind) const {
        return this->_elements[this->position(ind)]; }
      
      /*! \brief A reference to the element indexed by \a ind. */
      const R& operator() (const array<size_type> ind) {
        return this->_elements[this->position(ind)]; }
    
     protected:
      size_type position(const array<size_type>& ind) const {
        size_type result=ind[0];
        for(size_type d=1; d!=this->rank(); ++d) { result=this->_sizes[d-1]+ind[d]; }
        return result;
      }
      static array<size_type> compute_strides(const array<size_type>& sz) {
        size_type rnk=sz.size(); array<size_type> result(rnk); result[rnk-1]=1;
        for(size_type d=rnk; d!=0; --d) { result[d-1]=sz[d]*result[d]; }
        return result;
      }
     private:
      array<size_type> _sizes;
      array<size_type> _strides;
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
