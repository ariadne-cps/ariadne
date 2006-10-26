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
#include "../linear_algebra/multi_index.h"


namespace Ariadne {
  namespace LinearAlgebra {
    
      
    /*! \ingroup LinearAlgebra
     *  \brief Base class for tensor expressions. 
     */
    template<class E> 
    class TensorExpression {
     public:
      typedef IndexArray index_array_type;
      
//      typename E::scalar_type operator()(const index_array_type& i) const {
//        return static_cast<const E&>(*this)(i);
//      }
    };
    
    /*! \ingroup LinearAlgebra
     *  \brief Base class for symmetric tensor expressions. 
     */
    template<class E> 
    class SymmetricTensorExpression 
      : public TensorExpression<E>
    {
      typedef MultiIndex multi_index_type;
      
      const E& operator() () const { return static_cast<const E&>(*this); }
      E& operator() () { return static_cast<E&>(*this); }
        
//      typename E::scalar_type operator()(const size_type& i, const multi_index_type& j) const {
//        return static_cast<const E&>(*this)(i,j);
//      }
    };
    
    /*! \ingroup LinearAlgebra
     *  \brief Base class for derivative tensor expressions. 
     */
    template<class E> 
    class DerivativeTensorExpression 
      : public TensorExpression<E>
    {
      typedef MultiIndex multi_index_type;
      
      const E& operator() () const { return static_cast<const E&>(*this); }
      E& operator() () { return static_cast<E&>(*this); }
        
//      typename E::scalar_type operator()(const size_type& i, const multi_index_type& j) const {
//        return static_cast<const E&>(*this)(i,j);
//      }
    };
    
  
    /*! \ingroup LinearAlgebra
     *  \brief A tensor. 
     */
    template<class R> 
    class Tensor {
     public:
      typedef array<size_type> index_array_type;
      typedef R scalar_type;
     
      /*! \brief Construct a rank-three tensor with sizes \a m, \a n1 and \a n2. */
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
      
      /*! \brief A constant reference to the element indexed by \a ind. */
      const R& operator() (const index_array_type ind) const {
        return this->_elements[this->position(ind)]; }
      
      /*! \brief A reference to the element indexed by \a ind. */
      R& operator() (const index_array_type ind) {
        return this->_elements[this->position(ind)]; }
    
      //@{
      //! \name Element acces for low-rank tensors
      /*! \brief A constant reference to the \a i-th element of a rank-1 tensor. */
      const R& operator() (const size_type& i) const {
        assert(this->rank()==1);
        return this->_elements[i];
      }
      
      /*! \brief A reference to the \a i-th element of a rank-1 tensor. */
      R& operator() (const size_type& i) {
        assert(this->rank()==1);
        return this->_elements[i];
      }
      
      /*! \brief A constant reference to the (\a i,\a j)-th element of a rank-2 tensor. */
      const R& operator() (const size_type& i, const size_type& j) const {
        assert(this->rank()==2);
        return this->_elements[i*this->_sizes[1]+j];
      }
      
      /*! \brief A reference to the (\a i,\a j)-th element of a rank-2 tensor. */
      R& operator() (const size_type& i, const size_type& j) {
        assert(this->rank()==2);
        return this->_elements[i*this->_sizes[1]+j];
      }
      
      /*! \brief A constant reference to the (\a i,\a j,\a k)-th element of a rank-3 tensor. */
      const R& operator() (const size_type& i, const size_type& j, const size_type& k) const {
        assert(this->rank()==3);
        return this->_elements[(i*_sizes[1]+j)*_sizes[2]+k];
      }
      
      /*! \brief A reference to the (\a i,\a j,\a k)-th element element of a rank-3 tensor. */
      R& operator() (const size_type& i, const size_type& j, const size_type& k) {
        assert(this->rank()==3);
        return this->_elements[(i*_sizes[1]+j)*_sizes[2]+k];
      }
      //@}
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
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
  
    
    
    /*!\ingroup LinearAlgebra
     * \brief A tensor representing a derivative, with a single covariant (upper)
     * index, and symmetric contravariant (lower) indices. 
     */
    template<class R> 
    class SymmetricTensor 
      : public SymmetricTensorExpression< SymmetricTensor<R> >
    {
     public:
      typedef R scalar_type;
     
      typedef MultiIndex multi_index_type;
      typedef IndexArray index_array_type;
    
      /*! \brief Construct a tensor with 1 upper index of size m, and deg lower indices of size n. */
      SymmetricTensor(size_type n, size_type deg)
        : _arg_size(n), _degree(deg), _elements(choose(n+deg-1,deg))
      { 
        //std::cerr << "SymmetricTensor(size_type m, size_type n, size_type deg)" << std::endl;
      }
      
      /*! \brief Construct a tensor with 1 upper index of size m, and deg lower indices of size n. */
      SymmetricTensor(size_type n, size_type deg, const R* ptr) 
        : _arg_size(n), _degree(deg), _elements(choose(n+deg-1,deg))
      { _elements.fill(ptr); }
      
      /*! \brief Construct from a string. */
      SymmetricTensor(const std::string& str);
      
      /*! \brief Copy constructor. */
      template<class E>
      SymmetricTensor(const SymmetricTensorExpression<E>& t)
        : _arg_size(t().space_size()), _degree(t().degree()), 
          _elements(t().number_of_elements()) 
      { 
          MultiIndexIterator curr(this->_arg_size,this->_degree);
          MultiIndexIterator end(this->_arg_size,this->_degree+1);
          for( ; curr!=end; ++curr) {
            this->_elements[curr->position()]=t(*curr);
          }
      }

      /*! \brief Copy assignment operator. */
      template<class E>
      SymmetricTensor<R>& operator=(const SymmetricTensorExpression<E>& t) {
        if(this!=&t) {
          const E& e=t();
          this->_arg_size=e.argument_size();
          this->_degree=e.degree();
          this->_elements().resize(this->number_of_elements());
          MultiIndexIterator curr(this->_arg_size,this->_degree);
          MultiIndexIterator end(this->_arg_size,this->_degree+1);
          for( ; curr!=end; ++curr) {
            this->_elements[curr->position()]=t(*curr);
          }
        }     
        return *this;
      }        
      
      /*! \brief The rank of the tensor. */
      size_type rank() const { return this->degree(); }
      
      /*! \brief The number of lower indices. */
      size_type degree() const { return this->_degree; }
      
      /*! \brief The rank of the tensor. */
      size_type number_of_elements() const { 
        return this->_arg_size^this->_degree; }
        
      /*! \brief The size of the lower indices. */
      const size_type& argument_size() const { return this->_arg_size; }

      /*! \brief The size of the \a i th dimension. */
      size_type size(const size_type& i) const { return this->_arg_size; }
      
      /*! \brief A constant reference to the element indexed by \a j. */
      const R& operator() (const multi_index_type j) const {
        return this->_elements[this->position(j)]; }
      
      /*! \brief A reference to the element indexed by \a j. */
      R& operator() (const multi_index_type j) {
        //std::cerr << "SymmetricTensor<R>::operator() (const multi_index_type j)" << std::endl;
        return this->_elements[this->position(j)]; }
    
      /*! \brief A constant reference to the element indexed by \a i. */
      const R& operator() (const index_array_type i) const {
        return this->_elements[this->position(i)]; }
      
      /*! \brief A reference to the element indexed by \a i. */
      R& operator() (const index_array_type i) {
        return this->_elements[this->position(i)]; }
    
       /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     protected:
      size_type position(const multi_index_type& j) const;
      size_type position(const index_array_type& i) const;
      static size_type number_of_elements(const size_type& n, const size_type& d) {
        return choose(n+d-1,d); }
     private:
      size_type _arg_size;
      size_type _degree;
      array<R> _elements;
    };
  
    template<class R>
    size_type 
    SymmetricTensor<R>::position(const multi_index_type& j) const
    {
      assert(j.number_of_variables()==this->_arg_size);
      assert(j.degree()==this->_degree);
      return j.position();
    }
    
    template<class R>
    size_type 
    SymmetricTensor<R>::position(const index_array_type& i) const
    {
      assert(i.size()==this->_degree);
      return multi_index_type(this->_arg_size,i).position();
    }
    



    /*!\ingroup LinearAlgebra
     * \brief A tensor representing a derivative, with a single covariant (upper)
     * index, and symmetric contravariant (lower) indices. 
     */
    template<class R> 
    class DerivativeTensor 
      : public DerivativeTensorExpression< DerivativeTensor<R> >
    {
     public:
      typedef R scalar_type;
     
      typedef MultiIndex multi_index_type;
      typedef IndexArray index_array_type;
    
      /*! \brief Construct a tensor with 1 upper index of size m, and deg lower indices of size n. */
      DerivativeTensor(size_type m, size_type n, size_type deg)
        : _res_size(m), _arg_size(n), _degree(deg), _elements(m*choose(n+deg-1,deg))
      { 
        //std::cerr << "DerivativeTensor(size_type m, size_type n, size_type deg)" << std::endl;
      }
      
      /*! \brief Construct a tensor with 1 upper index of size m, and deg lower indices of size n. */
      DerivativeTensor(size_type m, size_type n, size_type deg, const R* ptr) 
        : _res_size(m), _arg_size(n), _degree(deg), _elements(m*choose(n+deg-1,deg))
      { _elements.fill(ptr); }
      
      /*! \brief Construct from a string. */
      DerivativeTensor(const std::string& str);
      
      /*! \brief Construct from a vector. */
      DerivativeTensor(const Vector<R>& v);
      
      /*! \brief Construct from a matrix. */
      DerivativeTensor(const Matrix<R>& v);
      
      /*! \brief Copy constructor. */
      template<class E>
      DerivativeTensor(const DerivativeTensorExpression<E>& t)
        : _res_size(t().result_size()), _arg_size(t().argument_size()), _degree(t().degree()), 
          _elements(t().number_of_elements()) 
      { 
          // TODO: Assign elements        
      }

      /*! \brief Copy assignment operator. */
      template<class E>
      DerivativeTensor<R>& operator=(const DerivativeTensorExpression<E>& t) {
        if(this!=&t) {
          const E& e=t();
          this->_res_size=e.result_size();
          this->_arg_size=e.argument_size();
          this->_degree=e.degree();
          this->_elements().resize(this->number_of_elements());
          MultiIndexIterator curr(this->_res_size,this->_arg_size,this->_degree);
          MultiIndexIterator end(this->_res_size,this->_arg_size,this->_degree+1);
          for( ; curr!=end; ++curr) {
            this->_elements[curr->position()]=t(*curr);
          }
        }     
        return *this;
      }        
      
      /*! \brief The rank of the tensor. */
      size_type rank() const { return this->degree()+1; }
      
      /*! \brief The number of lower indices. */
      size_type degree() const { return this->_degree; }
      
      /*! \brief The rank of the tensor. */
      size_type number_of_elements() const { 
        return this->_res_size*this->_arg_size^this->_degree; }
        
      /*! \brief The size of the upper index. */
      const size_type& result_size() const { return this->_res_size; }
      /*! \brief The size of the lower indices. */
      const size_type& argument_size() const { return this->_arg_size; }

      /*! \brief The size of the \a i th dimension. */
      size_type size(const size_type& i) const { return i==0 ? this->_res_size : this->_arg_size; }
      
      /*! \brief A constant reference to the element indexed by \a (i,j). */
      const R& operator() (const size_type& i, const multi_index_type j) const {
        assert(j.degree()==this->degree());
        return this->_elements[this->position(i,j)]; }
      
      /*! \brief A reference to the element indexed by \a (i,j). */
      R& operator() (const size_type& i, const multi_index_type j) {
        //std::cerr << "DerivativeTensor<R>::operator() (const size_type& i, const multi_index_type j)" << std::endl;
        assert(j.degree()==this->degree());
        return this->_elements[this->position(i,j)]; }
    
      /*! \brief A constant reference to the element indexed by \a i. */
      //const R& operator() (const size_type& i, const index_array_type i) const {
      //  return this->_elements[this->position(i)]; }
      
      /*! \brief A reference to the element indexed by \a ind. */
      //R& operator() (const index_array_type i) {
        //return this->_elements[this->position(i)]; }
    
      static DerivativeTensor<R> product(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2);
      static DerivativeTensor<R> product(const DerivativeTensor<R>& T, const Matrix<R>& A);
      static DerivativeTensor<R> product(const DerivativeTensor<R>& T, const Vector<R>& v);
        
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     protected:
      size_type position(const size_type& i, const multi_index_type& j) const;
      size_type position(const size_type& i, const index_array_type& i) const;
      static size_type number_of_elements(const size_type& nr, const size_type& na, const size_type& d) {
        return nr*choose(na+d-1,d); }
     private:
      size_type _res_size;
      size_type _arg_size;
      size_type _degree;
      array<R> _elements;
    };
  
    template<class R>
    size_type 
    DerivativeTensor<R>::position(const size_type& i, const multi_index_type& j) const
    {
      assert(i<this->_res_size);
      assert(j.number_of_variables()==this->_arg_size);
      assert(j.degree()==this->_degree);
      return i*choose(this->_arg_size+this->_degree-1,this->_degree)+j.position();
    }
    
    template<class R>
    DerivativeTensor<R>
    operator*(const DerivativeTensor<R>& T, const Vector<R>& v) {
      return DerivativeTensor<R>::product(T,v);
    }
    
    template<class R>
    DerivativeTensor<R>
    operator*(const DerivativeTensor<R>& T, const Matrix<R>& A) {
      return DerivativeTensor<R>::product(T,A);
    }

    template<class R>
    DerivativeTensor<R>
    operator*(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2) {
      return DerivativeTensor<R>::product(T1,T2);
    }

    
    
    template<class R>
    inline
    std::ostream&
    operator<< (std::ostream& os, const Tensor<R>& T) {
      return T.write(os);
    } 

    template<class R>
    inline
    std::ostream&
    operator<< (std::ostream& os, const SymmetricTensor<R>& T) {
      return T.write(os);
    } 
    
    template<class R>
    inline
    std::ostream&
    operator<< (std::ostream& os, const DerivativeTensor<R>& T) {
      return T.write(os);
    } 
    
  }
}

#endif /* _ARIADNE_TENSOR_H */
