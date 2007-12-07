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

#ifndef ARIADNE_TENSOR_H
#define ARIADNE_TENSOR_H

#include <iosfwd>

#include "../base/array.h"
#include "../numeric/traits.h"
#include "../linear_algebra/exceptions.h"
#include "../linear_algebra/multi_index.h"
#include "../linear_algebra/tensor_expression.h"


namespace Ariadne {
  namespace LinearAlgebra {
    
    // Forward declarations
    template<class R> class Vector;
    template<class R> class Matrix;
      
  
    /*! \ingroup LinearAlgebra
     *  \brief A tensor. 
     */
    template<class R> 
    class Tensor 
      : public TensorExpression< Tensor<R> >
    {
     public:
      typedef array<size_type> index_array_type;
      typedef R scalar_type;
     
      /*! \brief Construct a tensor with sizes \a sz. */
      Tensor(const array<size_type>& sz);
      
      /*! \brief Construct a tensor with sizes \a sz from the elements beginning at \a ptr. */
      Tensor(const array<size_type>& sz, const R* ptr);
      
      /*! \brief Conversion from a tensor over another real type. */
      template<class Rl> Tensor(const Tensor<Rl>& t);
      
      /*! \brief Copy constructor. */
      Tensor(const Tensor<R>& t);

      /*! \brief Copy assignment operator. */
      Tensor<R>& operator=(const Tensor<R>& t);
      
      /*! \brief The rank of the tensor. */
      size_type rank() const;
      
      /*! \brief The rank of the tensor. */
      size_type number_of_elements() const;
        
      /*! \brief The array of sizes in each dimension. */
      const array<size_type>& sizes() const;
      /*! \brief The size of the \a i th dimension. */
      size_type size(const size_type& i) const;
      
      /*! \brief A constant reference to the element indexed by \a ind. */
      const R& operator() (const index_array_type ind) const;
      
      /*! \brief A reference to the element indexed by \a ind. */
      R& operator() (const index_array_type ind);
    
      /*! \brief A pointer to the beginning of the array of elements. */
      R* begin();
      
      /*! \brief A constant pointer to the end of the array of elements. */
      const R* begin() const;
      
      //@{
      //! \name Constructors for low-rank tensors
      /*! \brief Construct a rank-three tensor with sizes \a m, \a n1 and \a n2. */
      Tensor(size_type m, size_type n1, size_type n2);
      //@}
        
      //@{
      //! \name Element acces for low-rank tensors
      /*! \brief A constant reference to the \a i-th element of a rank-1 tensor. */
      const R& operator() (const size_type& i) const;
      
      /*! \brief A reference to the \a i-th element of a rank-1 tensor. */
      R& operator() (const size_type& i);
      
      /*! \brief A constant reference to the (\a i,\a j)-th element of a rank-2 tensor. */
      const R& operator() (const size_type& i, const size_type& j) const;
      
      /*! \brief A reference to the (\a i,\a j)-th element of a rank-2 tensor. */
      R& operator() (const size_type& i, const size_type& j);
      
      /*! \brief A constant reference to the (\a i,\a j,\a k)-th element of a rank-3 tensor. */
      const R& operator() (const size_type& i, const size_type& j, const size_type& k) const;
      
      /*! \brief A reference to the (\a i,\a j,\a k)-th element element of a rank-3 tensor. */
      R& operator() (const size_type& i, const size_type& j, const size_type& k);
      //@}
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     protected:
      size_type position(const array<size_type>& ind) const;
      static array<size_type> compute_strides(const array<size_type>& sz);
     private:
      static void _instantiate();
     private:
      array<size_type> _sizes;
      array<size_type> _strides;
      array<R> _elements;
    };
  
    template<class R>
    Tensor<typename Numeric::traits<R>::arithmetic_type> 
    operator+(const Tensor<R>& T1, const Tensor<R>& T2);
    
    template<class R>
    Tensor<typename Numeric::traits<R>::arithmetic_type> 
    operator*(const Tensor<R>& T, const R& s);
    
    template<class R>
    std::ostream& operator<<(std::ostream& os, const Tensor<R>& T);
    
    
    
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
      SymmetricTensor(size_type n, size_type deg);
      
      /*! \brief Construct a tensor with 1 upper index of size m, and deg lower indices of size n. */
      SymmetricTensor(size_type n, size_type deg, const R* ptr);
      
      /*! \brief Copy constructor. */
      template<class E> SymmetricTensor(const SymmetricTensorExpression<E>& t);

      /*! \brief Copy assignment operator. */
      template<class E> SymmetricTensor<R>& operator=(const SymmetricTensorExpression<E>& t);
    
      /*! \brief The rank of the tensor. */
      size_type rank() const;
      
      /*! \brief The number of lower indices. */
      size_type degree() const;
      
      /*! \brief The total number of elements (note that many are the same). */
      size_type number_of_elements() const;
        
      /*! \brief The total number of elements (note that many are the same). */
      size_type number_of_independent_elements() const;
        
      /*! \brief The size of the lower indices. */
      const size_type& argument_size() const;

      /*! \brief The size of the \a i th dimension (must be the same in all dimensions). */
      size_type size(const size_type& i) const;
      
      /*! \brief A constant reference to the element indexed by \a j. */
      const R& operator() (const multi_index_type j) const;
      
      /*! \brief A reference to the element indexed by \a j. */
      R& operator() (const multi_index_type j);
    
      /*! \brief A constant reference to the element indexed by \a i. */
      const R& operator() (const index_array_type i) const;
      
      /*! \brief A reference to the element indexed by \a i. */
      R& operator() (const index_array_type i);
      
      /*! \brief A pointer to the beginning of the array of elements. */
      R* begin();
      
      /*! \brief A constant pointer to the end of the array of elements. */
      const R* begin() const;
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     protected:
      size_type position(const multi_index_type& j) const;
      size_type position(const index_array_type& i) const;
      static size_type number_of_independent_elements(const size_type& n, const size_type& d);
     private:
      static void _instantiate();
     private:
      size_type _arg_size;
      size_type _degree;
      array<R> _elements;
    };

    template<class R>
    SymmetricTensor<typename Numeric::traits<R>::arithmetic_type> 
    operator+(const SymmetricTensor<R>& T1, const SymmetricTensor<R>& T2);
    
    template<class R>
    SymmetricTensor<typename Numeric::traits<R>::arithmetic_type> 
    operator*(const SymmetricTensor<R>& T, const R& s);
    
    template<class R> std::ostream& operator<<(std::ostream& os, const SymmetricTensor<R>& S); 
    


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
      DerivativeTensor(size_type m, size_type n, size_type deg);
      
      /*! \brief Construct a tensor with 1 upper index of size m, and deg lower indices of size n. */
      DerivativeTensor(size_type m, size_type n, size_type deg, const R* ptr);
      
      /*! \brief Construct from a vector. */
      DerivativeTensor(const Vector<R>& v);
      
      /*! \brief Construct from a matrix. */
      DerivativeTensor(const Matrix<R>& v);
      
      /*! \brief Copy constructor. */
      template<class E> DerivativeTensor(const DerivativeTensorExpression<E>& t);
      
      /*! \brief Copy assignment operator. */
      template<class E> DerivativeTensor<R>& operator=(const DerivativeTensorExpression<E>& t);
      
      /*! \brief The rank of the tensor. */
      size_type rank() const;
      
      /*! \brief The number of lower indices. */
      size_type degree() const;
      
      /*! \brief The total number of elements. */
      size_type number_of_elements() const;
        
      /*! \brief The number of independent scalar values needed to define the tensor. */
      size_type number_of_independent_elements() const;
        
      /*! \brief The size of the upper index. */
      const size_type& result_size() const;
      /*! \brief The size of the lower indices. */
      const size_type& argument_size() const;

      /*! \brief The size of the \a i th dimension. */
      size_type size(const size_type& i) const;
      
      /*! \brief A constant reference to the element indexed by \a (i,j). */
      const R& operator() (const size_type& i, const multi_index_type j) const;
      
      /*! \brief A reference to the element indexed by \a (i,j). */
      R& operator() (const size_type& i, const multi_index_type j);
    
      /*! \brief A constant reference to the element indexed by \a i. */
      //const R& operator() (const size_type& i, const index_array_type i) const;
      
      /*! \brief A reference to the element indexed by \a ind. */
      //R& operator() (const index_array_type i);
    
      /*! \brief A pointer to the beginning of the array of elements. */
      R* begin();
      
      /*! \brief A constant pointer to the end of the array of elements. */
      const R* begin() const;
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     protected:
      template<class E> void assign(const E& e);
      size_type position(const size_type& i, const multi_index_type& j) const;
      size_type position(const size_type& i, const index_array_type& i) const;
      static size_type number_of_independent_elements(const size_type& nr, const size_type& na, const size_type& d);
     private:
      static void _instantiate();
     private:
      size_type _res_size;
      size_type _arg_size;
      size_type _degree;
      array<R> _elements;
    };
  
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator+(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2);
    
    template<class R> 
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator-(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2);
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T, const R& s);
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T, const Vector<R>& v);
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T, const Matrix<R>& A);
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2);
    
    template<class R> std::ostream& operator<<(std::ostream& os, const DerivativeTensor<R>& D); 

  }
}

#include "tensor.inline.h"

#endif /* ARIADNE_TENSOR_H */
