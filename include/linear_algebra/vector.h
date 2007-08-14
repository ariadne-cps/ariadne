/***************************************************************************
 *            vector.h
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
 
/*! \file vector.h
 *  \brief Vector types and vector operations.
  */

#ifndef ARIADNE_VECTOR_H
#define ARIADNE_VECTOR_H 

#include <iosfwd>
#include <algorithm>

#include "../base/types.h"
#include "../base/array.h"
#include "../numeric/integer.h"
#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/exceptions.h"
#include "../linear_algebra/vector_expression.h"

namespace Ariadne {
  namespace LinearAlgebra {
    
    

    /*! \ingroup LinearAlgebra
     *  \brief A vector over \a R. 
     */
    template<class R>
    class Vector
      : public VectorExpression< Vector<R> >
    {
     private:
      array<R> _array;
     public:
      typedef R value_type;

      /*! \brief Construct a vector of size 0. */
      explicit Vector();
      /*! \brief Construct the zero vector of size \a n. */
      explicit Vector(const size_type& n);
      /*! \brief Construct the zero vector of size \a n with all elements initialized to \a x. */
      explicit Vector(const size_type& n, const R& x);
      /*! \brief Construct a vector from the array \a ary. */
      explicit Vector(const array<R>& ary);
      /*! \brief Construct a vector of size \a n from the array beginning at \a ptr. */
      explicit Vector(const size_type& n, const R* ptr, const size_type& inc=1);

      /* Convert from a vector expression. */
      template<class E> Vector(const VectorExpression<E>& ve);
        
      /*! \brief Construct from a string literal of the form 
       *  "[v1,v2,...,vn]", where the vi are numerical literals for the vector elements. 
       */
      explicit Vector(const std::string& s); 

      /*! \brief Copy constructor. */
      Vector(const Vector<R>& v);
      /*! \brief Copy assignment operator. */
      Vector<R>& operator=(const Vector<R>& v);
      
      /*! \brief The equality operator. */
      bool operator==(const Vector<R>& v) const;
        
      /*! \brief The inequality operator. */
      bool operator!=(const Vector<R>& v) const;
      
      /*! \brief Returns the zero vector of size \a n. */
      static Vector<R> zero(const size_type& n);
      /*! \brief Returns the vector of size \a n whose elements are all 1. */
      static Vector<R> one(const size_type& n);
      /*! \brief Returns the unit vector of size \a n whose \a i th element is 1. */
      static Vector<R> unit(dimension_type n, dimension_type i);
        
      /*! \brief A reference to the array storing the elements data. */
      array<R>& data();
      /*! \brief The array storing the elements data. */
      const array<R>& data() const;

      /*! \brief Resize the vector to hold \a n elements. */
      void resize(const size_type& n);

      /*! \brief The number of elements of the vector. */
      size_type size() const;
      /*! \brief A pointer to the beginning of the data. */
      R* begin();
      /*! \brief A constant pointer to the beginning of the data. */
      const R* begin() const;
      /*! \brief A pointer to the end of the data. */
      R* end();
      /*! \brief A constant pointer to the end of the data. */
      const R* end() const;

      /*! \brief The increment between elements in the storage array. */
      size_type increment() const;

      /*! \brief A reference to the \a i th element. */
      R& operator() (const size_type& i);
      /*! \brief A constant reference to the \a i th element. */
      const R& operator() (const size_type& i) const;

      /*! \brief A reference to the \a i th element. */
      R& operator[] (const size_type& i);
      /*! \brief A constant reference to the \a i th element. */
      const R& operator[] (const size_type& i) const;

            
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

      /*! \brief The inner product of \a v1 and \a v2. */
      friend R inner_product<>(const Vector<R>& v1, const Vector<R>& v2);
      /*! \brief The supremum norm of \a v. */
      friend R LinearAlgebra::sup_norm<>(const Vector<R>& v);
      /*! \brief The supremum norm of \a v. */
      friend R LinearAlgebra::norm<>(const Vector<R>& v);
      /*! \brief The direct sum (concatentation) of v1 and v2. */
      friend Vector<R> direct_sum<R>(const Vector<R>& v1, const Vector<R>& v2);
      /*! \brief The concatentation of v1 and v2. */
      friend Vector<R> concatenate<R>(const Vector<R>& v1, const Vector<R>& v2);
      /*! \brief The concatentation of v and s. */
      friend Vector<R> concatenate<R>(const Vector<R>& v1, R& s);
#endif 

      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);
    };
    
    

    /*! \ingroup LinearAlgebra
     *  \brief A slice through an array or vector, with equally spaces increments. 
     */
    template<class R>
    class VectorSlice : public VectorExpression< VectorSlice<R> >
    {
     public:
      typedef R value_type;
     
      explicit VectorSlice(const size_type& size, R* begin, const size_type& increment=1u);
      VectorSlice(const Vector<R>& v);
      VectorSlice(const array<R>& a);
      
      size_type size() const;
      const R* begin() const;
      const R* end() const;
      size_type increment() const;
      
      const R& operator() (const size_type& i) const;
      R& operator() (const size_type& i);
     
      const R& operator[] (const size_type& i) const;
      R& operator[] (const size_type& i);
     
      template<class E> VectorSlice<R>& operator=(const VectorExpression< E >& v);
      
      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
     private:
      size_type _size;
      R* _begin;
      size_type _increment;
    };
    
    
    template<class R1, class E2> Vector<R1>& operator+=(Vector<R1>& v1, const VectorExpression<E2>& e2);
    template<class R1, class E2> Vector<R1>& operator-=(Vector<R1>& v1, const VectorExpression<E2>& e2);
    template<class R1, class R2> Vector<R1>& operator*=(Vector<R1>& v1, const R2& s2);
    template<class R1, class R2> Vector<R1>& operator/=(Vector<R1>& v1, const R2& s2);

    template<class E1, class E2> BinaryVectorVectorExpression<E1,E2,plus> operator+(const VectorExpression<E1>& e1, const VectorExpression<E2>& e2);
    template<class E1, class E2> BinaryVectorVectorExpression<E1,E2,minus> operator-(const VectorExpression<E1>& e1, const VectorExpression<E2>& e2);
    template<class E1, class E2> BinaryVectorScalarExpression<E1,E2,times> operator*(const E2& e2, const VectorExpression<E1>& e1);
    template<class E1, class E2> BinaryVectorScalarExpression<E1,E2,times> operator*(const VectorExpression<E1>& e1, const E2& e2);
    template<class E1, class E2> BinaryVectorScalarExpression<E1,typename Numeric::traits<E2>::closure_type,divides> operator/(const VectorExpression<E1>& e1, const E2& e2);



  template<class R> Vector<R> zero_vector(size_type n);
  template<class R> Vector<R> unit_vector(size_type n, size_type i);
	
	
  template<class R> Vector<R> midpoint(const Vector< Numeric::Interval<R> >& iv); 
  template<class R> bool encloses(const Vector< Numeric::Interval<R> >& iv,const Vector<R>& v); 
  template<class R> bool refines(const Vector< Numeric::Interval<R> >& iv1, const Vector< Numeric::Interval<R> >& iv2); 

  template<class R1, class R2> Vector<R1> approximation(const Vector<R2>& v);

  template<class R> std::ostream& operator<<(std::ostream& os, const Vector<R>& v);
  template<class R> std::istream& operator>>(std::istream& is, Vector<R>& v);
  template<class R> std::ostream& operator<<(std::ostream& os, const VectorSlice<R>& vs);
    
    
    
    
    
    template<class R> Vector<R> operator-(const Vector<R>& v);
    template<class R1, class R2> Vector<typename Numeric::traits<R1,R2>::arithmetic_type> operator+(const Vector<R1>& v1, const Vector<R2>& v2);
    template<class R1, class R2> Vector<class Numeric::traits<R1,R2>::arithmetic_type> operator-(const Vector<R1>& v1, const Vector<R2>& v2);
    template<class R1, class R2> Vector<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const R1& s, const Vector<R2>& v);
    template<class R1, class R2> Vector<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const Vector<R1>& v, const R2& s);
    template<class R1, class R2> Vector<typename Numeric::traits<R1,R2>::arithmetic_type> operator/(const Vector<R1>& v, const R2& s);
   
    
    template<class R> Vector<R> add_approx(const Vector<R>& u, const Vector<R>& v);
    template<class R> Vector<R> sub_approx(const Vector<R>& u, const Vector<R>& v);
    template<class R> Vector<R> mul_approx(const R& s, const Vector<R>& v);
    template<class R> Vector<R> mul_approx(const Vector<R>& v, const R& s);
    template<class R> Vector<R> div_approx(const Vector<R>& v, const R& s);
    
    
    template<class R> R inner_product(const Vector<R>& u, const Vector<R>& v);
    template<class R> Vector<R> direct_sum(const Vector<R>& v1, const Vector<R>& v2);
    template<class R> R sup_norm(const Vector<R>& v);
    template<class R> R norm(const Vector<R>& v);

    template<class R> Vector<R> concatenate(const Vector<R>& v1, const Vector<R>& v2);
    template<class R> Vector<R> concatenate(const Vector<R>& v, const R& s);
  }
}

#include "vector.inline.h"

#endif /* ARIADNE_VECTOR_H */
