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
 *  \brief vector_types of intervals.
  */

#ifndef _ARIADNE_INTERVAL_VECTOR_H
#define _ARIADNE_INTERVAL_VECTOR_H 

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../declarations.h"
#include "../numeric/interval.h"
#include "../linear_algebra/vector.h"

namespace Ariadne {
  namespace LinearAlgebra {
    
    /*! \brief A vector of intervals. */
    template<typename R>
    class interval_vector : public boost::numeric::ublas::vector< Interval<R> >
    {
     private:
      typedef boost::numeric::ublas::vector< Interval<R> > _Base;
     public:
      template<typename E> interval_vector(const boost::numeric::ublas::vector_expression<E> v) :_Base(v()) { }
      interval_vector() : _Base() { }
      interval_vector(const size_type& n) : _Base(n) { }
      interval_vector(const vector<R>& v);
      interval_vector(const vector<R>& v, const R& r);
      
      vector<R> centre() const;
      R radius() const;
      
      Interval<R> norm() const;
      R upper_norm() const;
    };
      
    template <typename R> interval_vector<R> operator+(const interval_vector<R>& iv1, const interval_vector<R> iv2);
    template <typename R> interval_vector<R> operator*(const Interval<R>& s, const interval_vector<R>& iv);
    template <typename R> interval_vector<R> operator*(const R& s, const interval_vector<R>& iv);

    
    template <typename R>
    inline
    interval_vector<R>::interval_vector(const vector<R>& v)
      : _Base(v.size()) 
    { 
      for(size_type i=0; i!=this->size(); ++i) {
        _Base::operator()(i) = Interval<R>(v(i),v(i));
      }
    }
    
    template <typename R>
    inline
    interval_vector<R>::interval_vector(const vector<R>& v, const R& r)
      : _Base(v.size()) 
    { 
      for(size_type i=0; i!=this->size(); ++i) {
        _Base::operator()(i) = Interval<R>(v(i)-r,v(i)+r);
      }
    }
      
    template <typename R>
    inline
    vector<R> 
    interval_vector<R>::centre() const
    {
      vector<R> result(this->size());
      const interval_vector<R>& v=*this;
      for(size_type i=0; i!=v.size(); ++i) {
        result(i) = v(i).centre();
      }
      return result;
    }
    
    template <typename R>
    inline
    R
    interval_vector<R>::radius() const
    {
      const interval_vector<R>& v=*this;
      R result=0;
      for(size_type i=0; i!=v.size(); ++i) {
        result = max(result,v(i).length());
      }
      return result;
    }
    
    template<typename R>
    inline
    Interval<R>
    interval_vector<R>::norm() const 
    {
      const interval_vector<R>& v=*this;
      R lower_bound=0;
      R upper_bound=0;
      for (size_type i=0; i<v.size(); i++) {
        if(!(v(i).lower()<=0 && v(i).upper()>=0)) {
          lower_bound=min(lower_bound,R(min(abs(v(i).lower()),abs(v(i).upper()))));
        }
        upper_bound=max(upper_bound,R(max(abs(v(i).lower()),abs(v(i).upper()))));
      }
      return Interval<R>(lower_bound,upper_bound);
    }
    
    template<typename R>
    inline
    R
    interval_vector<R>::upper_norm() const 
    {
      const interval_vector<R>& v=*this;
      R upper_bound=0;
      for (size_type i=0; i<v.size(); i++) {
        upper_bound=max(upper_bound,R(max(abs(v(i).lower()),abs(v(i).upper()))));
      }
      return upper_bound;
    }
    
    
    
    template <typename R>
    inline
    interval_vector<R> 
    operator+(const interval_vector<R>& iv, const vector<R> v)
    {
      interval_vector<R> result(v.size());
      for(size_type i=0; i!=result.size(); ++i) {
        result(i)=v(i)+iv(i);
      }
      return result;
    }

    template <typename R>
    inline
    interval_vector<R> 
    operator+(const vector<R>& v, const interval_vector<R> iv)
    {
      interval_vector<R> result(v.size());
      for(size_type i=0; i!=result.size(); ++i) {
        result(i)=v(i)+iv(i);
      }
      return result;
    }

    template <typename R>
    inline
    interval_vector<R> 
    operator+(const interval_vector<R>& iv1, const interval_vector<R> iv2)
    {
      interval_vector<R> result(iv1.size());
      for(size_type i=0; i!=result.size(); ++i) {
        result(i)=iv1(i)+iv2(i);
      }
      return result;
    }

    template<typename R>
    inline
    interval_vector<R>
    operator*(const Interval<R>& s, const interval_vector<R>& v)
    {
      interval_vector<R> result(v.size());
      for(size_type i=0; i!=result.size(); ++i) {
        result(i)=s*v(i);
      }
      return result;
    }

    template<typename R>
    inline
    interval_vector<R>
    operator*(const R& s, const interval_vector<R>& v)
    {
      return Interval<R>(s)*v;
    }

    template<typename R>
    inline
    interval_vector<R>
    operator*(const Interval<R>& s, const vector<R>& v)
    {
      return s*interval_vector<R>(v); 
    }


    
    template<typename R>
    inline
    Interval<R>
    norm(const interval_vector<R>& v)
    {
      return v.norm();
    }

    template<typename R>
    inline
    R
    upper_norm(const interval_vector<R>& v)
    {
      return v.upper_norm();
    }
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const interval_vector<R>& v);
    
  }
}  

#endif /* _ARIADNE_INTERVAL_VECTOR_H */
