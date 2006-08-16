/***************************************************************************
 *            interval_vector.h
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
 
/*! \file interval_vector.h
 *  \brief Vectors of intervals.
  */

#ifndef _ARIADNE_INTERVAL_VECTOR_H
#define _ARIADNE_INTERVAL_VECTOR_H 

#include <iosfwd>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../declarations.h"
#include "../numeric/interval.h"
#include "../linear_algebra/vector.h"

namespace Ariadne {
  namespace LinearAlgebra {
    
    /*! \brief A vector of intervals. */
    template<typename R>
    class IntervalVector : public Vector< Interval<R> >
    {
     private:
      typedef Vector< Interval<R> > _Base;
     public:
      /*! \brief Construct an interval vector of size zero. */
      explicit IntervalVector() : _Base() { }
      /*! \brief Construct an interval vector of size \a n, all of whose entries are the interval \f$[0,0]\f$. */
      explicit IntervalVector(const size_type& n) : _Base(n) { 
        for(size_type i=0; i!=n; ++i) { (*this)(i)=Interval<R>(0); } }
      /*! \brief Construct an interval vector of size \a n from the array starting at \a ptr. */
      explicit IntervalVector(const size_type& n, const Interval<R>* ptr) : _Base(n,ptr) { }
      /*! \brief Construct an interval vector centred at \a v of radius \a r. */
      explicit IntervalVector(const Vector<R>& v, const R& r) : _Base(v.size()) 
      { 
        for(size_type i=0; i!=this->size(); ++i) {
          this->_Base::operator()(i) = v(i) + Interval<R>(-r,r); }
      }
      /*! \brief Convert from a vector. */
      IntervalVector(const Vector<R>& v) : _Base(v) { }
      
      /* Convert from a vector expression. */
      template<typename E> IntervalVector(const boost::numeric::ublas::vector_expression<E>& v) : _Base(v) { }
 
      /*! \brief Construct from a string literal of the form "[[l1,u1],[l2,u2],...,[ln,un]]". */
      explicit IntervalVector(const std::string& str) : _Base(str) { }

#ifdef DOXYGEN
      /*! \brief Copy constructor. */
      IntervalVector(const IntervalVector<R>& v);
      /*! \brief Copy assignment operator. */
      IntervalVector<R>& operator=(const IntervalVector<R>& v);
#endif
      /*! \brief The centre of the interval vector. */
      Vector<R> centre() const;
      /*! \brief The maximum distance of an element of the interval vector from the centre in the supremum norm. */
      R radius() const;
      
      /*! \brief The highest possible supremum norm. */
      R upper_norm() const;
    };
      
    template <typename R>
    inline
    Vector<R> 
    IntervalVector<R>::centre() const
    {
      Vector<R> result(this->size());
      const IntervalVector<R>& v=*this;
      for(size_type i=0; i!=v.size(); ++i) {
        result(i) = v(i).centre();
      }
      return result;
    }
    
    template <typename R>
    inline
    R
    IntervalVector<R>::radius() const
    {
      const IntervalVector<R>& v=*this;
      R result=0;
      for(size_type i=0; i!=v.size(); ++i) {
        result = max(result,v(i).length());
      }
      return result;
    }
    
    template<typename R>
    inline
    R
    IntervalVector<R>::upper_norm() const 
    {
      const IntervalVector<R>& v=*this;
      R upper_bound=0;
      for (size_type i=0; i<v.size(); i++) {
        upper_bound=max(upper_bound,abs(v(i)).upper());
      }
      return upper_bound;
    }
    
    template <typename R>
    inline
    IntervalVector<R> 
    operator+(const IntervalVector<R>& iv, const Vector<R>& v)
    {
      return iv+IntervalVector<R>(v);
    }

    template <typename R>
    inline
    IntervalVector<R> 
    operator+(const Vector<R>& v, const IntervalVector<R>& iv)
    {
      return IntervalVector<R>(v)+iv;
    }

    template <typename R>
    inline
    IntervalVector<R> 
    operator-(const IntervalVector<R>& iv, const Vector<R>& v)
    {
      return iv-IntervalVector<R>(v);
    }

    template <typename R>
    inline
    IntervalVector<R> 
    operator-(const Vector<R>& v, const IntervalVector<R>& iv)
    {
      return IntervalVector<R>(v)-iv;
    }

    template<typename R>
    inline
    IntervalVector<R>
    operator*(const R& s, const IntervalVector<R>& v)
    {
      return Interval<R>(s)*v;
    }

    template<typename R>
    inline
    IntervalVector<R>
    operator*(const Interval<R>& s, const Vector<R>& v)
    {
      return s*IntervalVector<R>(v); 
    }

    template<typename R>
    inline
    IntervalVector<R>
    operator*(const IntervalVector<R>& v, const R& s)
    {
      return Interval<R>(s)*v;
    }

    template<typename R>
    inline
    IntervalVector<R>
    operator*(const Vector<R>& v, const Interval<R>& s)
    {
      return s*IntervalVector<R>(v); 
    }

    template<typename R>
    inline
    IntervalVector<R>
    operator/(const IntervalVector<R>& v, const R& s)
    {
      return v/Interval<R>(s); 
    }

    template<typename R>
    inline
    IntervalVector<R>
    operator/(const Vector<R>& v, const Interval<R>& s)
    {
      return IntervalVector<R>(v)/s; 
    }

    
    template<typename R>
    inline
    Interval<R>
    norm(const IntervalVector<R>& v)
    {
      return v.norm();
    }

    template<typename R>
    inline
    R
    upper_norm(const IntervalVector<R>& v)
    {
      return v.upper_norm();
    }
    
    template <typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const IntervalVector<R>& v)
    {
      return v.write(os); 
    }
    
    template <typename R>
    inline
    std::istream&
    operator>>(std::istream& is, IntervalVector<R>& v)
    {
      return v.read(is); 
    }
    
    
  }
}  

#endif /* _ARIADNE_INTERVAL_VECTOR_H */
