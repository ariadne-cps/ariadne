/***************************************************************************
 *            covector.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file covector.h
 *  \brief Covectors.
  */

#ifndef ARIADNE_COVECTOR_H
#define ARIADNE_COVECTOR_H 

#include <iosfwd>
#include <algorithm>

#include "base/types.h"
#include "base/array.h"
#include "numeric/integer.h"
#include "numeric/float.h"
#include "numeric/interval.h"

#include "linear_algebra/exceptions.h"

namespace Ariadne {
  
    
    template<class R> class Vector;
    template<class R> class Covector;
    template<class R> class Matrix;

    /*! \ingroup LinearAlgebra
     *  \brief A covector over \a R. A covector is a dual to a vector, and is represented as a row vector.
     */
    template<class R>
    class Covector
    {
     private:
      array<R> _data;
     public:
      typedef R value_type;

      /*! \brief Construct a vector of size 0. */
      explicit Covector();
      /*! \brief Construct the zero vector of size \a n. */
      explicit Covector(const size_type& n);
      /*! \brief Construct the zero vector of size \a n with all elements initialized to \a x. */
      explicit Covector(const size_type& n, const R& x);
      /*! \brief Construct a vector from the array \a ary. */
      template<class RR> explicit Covector(const array<RR>& ary);
      /*! \brief Construct a vector of size \a n from the array beginning at \a ptr with increment inc. */
      template<class RR> explicit Covector(const size_type& n, const RR* ptr, const size_type& inc=1);

      /*! \brief Copy constructor. */
      Covector(const Covector<R>& v);
      /*! \brief Conversion constructor. */
      template<class RR> Covector(const Covector<RR>& v);
      /*! \brief Copy assignment operator. */
      Covector<R>& operator=(const Covector<R>& v);
      /*! \brief Copy assignment operator. */
      template<class RR> Covector<R>& operator=(const Covector<RR>& v);
      
      /*! \brief The equality operator. */
      bool operator==(const Covector<R>& v) const;
        
      /*! \brief The inequality operator. */
      bool operator!=(const Covector<R>& v) const;
      
      /*! \brief Returns the zero vector of size \a n. */
      static Covector<R> zero(const size_type& n);
      /*! \brief Returns the vector of size \a n whose elements are all 1. */
      static Covector<R> one(const size_type& n);
      /*! \brief Returns the unit vector of size \a n whose \a i th element is 1. */
      static Covector<R> unit(dimension_type n, dimension_type i);
        
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

     private:
      static void instantiate();
    };
    
    

    
    template<class R1, class R2> Covector<R1>& operator+=(Covector<R1>& v1, const Covector<R2>& e2);
    template<class R1, class R2> Covector<R1>& operator-=(Covector<R1>& v1, const Covector<R2>& e2);
    template<class R1, class R2> Covector<R1>& operator*=(Covector<R1>& v1, const R2& s2);
    template<class R1, class R2> Covector<R1>& operator/=(Covector<R1>& v1, const R2& s2);

    template<class R> Covector<R> operator-(const Covector<R>& v);
    template<class R1, class R2> Covector<typename traits<R1,R2>::arithmetic_type> operator+(const Covector<R1>& v1, const Covector<R2>& v2);
    template<class R1, class R2> Covector<class traits<R1,R2>::arithmetic_type> operator-(const Covector<R1>& v1, const Covector<R2>& v2);
    template<class R1, class R2> Covector<typename traits<R1,R2>::arithmetic_type> operator*(const R1& s, const Covector<R2>& v);
    template<class R1, class R2> Covector<typename traits<R1,R2>::arithmetic_type> operator*(const Covector<R1>& v, const R2& s);
    template<class R1, class R2> Covector<typename traits<R1,R2>::arithmetic_type> operator/(const Covector<R1>& v, const R2& s);
    template<class R1, class R2> typename traits<R1,R2>::arithmetic_type operator*(const Covector<R1>& cv, const Vector<R2>& v);
    template<class R1, class R2> Covector<typename traits<R1,R2>::arithmetic_type> operator*(const Covector<R1>& cv, const Matrix<R2>& A);
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const Vector<R1>& v, const Covector<R2>& cv);
   
    template<class R> std::ostream& operator<<(std::ostream& os, const Covector<R>& cv);


    template<class R> Covector<R> zero_covector(size_type n);
    template<class R> Covector<R> unit_covector(size_type n, size_type i);
	
	
    template<class R> Covector<R> midpoint(const Covector< Interval<R> >& iv); 
    template<class R> bool encloses(const Covector< Interval<R> >& iv,const Covector<R>& v); 
    template<class R> bool refines(const Covector< Interval<R> >& iv1, const Covector< Interval<R> >& iv2); 

    template<class R1, class R2> Covector<R1> approximation(const Covector<R2>& v);

    template<class R> std::ostream& operator<<(std::ostream& os, const Covector<R>& v);
    

    template<class R> Covector<R> add_approx(const Covector<R>& u, const Covector<R>& v);
    template<class R> Covector<R> sub_approx(const Covector<R>& u, const Covector<R>& v);
    template<class R> Covector<R> mul_approx(const R& s, const Covector<R>& v);
    template<class R> Covector<R> mul_approx(const Covector<R>& v, const R& s);
    template<class R> Covector<R> div_approx(const Covector<R>& v, const R& s);
    
    
    template<class R> Covector<R> direct_sum(const Covector<R>& v1, const Covector<R>& v2);
    template<class R> R sup_norm(const Covector<R>& v);
    template<class R> R norm(const Covector<R>& v);

    template<class R> Covector<R> concatenate(const Covector<R>& v1, const Covector<R>& v2);
    template<class R> Covector<R> concatenate(const Covector<R>& v, const R& s);
  
} // namespace Ariadne

#include "covector.inline.h"

#endif /* ARIADNE_COVECTOR_H */
