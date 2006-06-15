/***************************************************************************
 *            interval.h
 *
 *  Wed 4 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file interval.h
 *  \brief Intervals of real number types (currently implemented using Boost).
 */
 
#ifndef _ARIADNE_INTERVAL_H
#define _ARIADNE_INTERVAL_H

#include <iostream>

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

#include "../declarations.h"

/* No input routine for intervals defined by boost */
namespace boost {
  namespace numeric {
    template<typename R>
    std::istream&
    operator>>(std::istream& is, interval<R>& ivl)
    {
      R l,u;
      char c1,c2,c3;
      is >> c1 >> l >> c2 >> u >> c3;
      ivl=interval<R>(l,u);
      return is;
    }
  }
}


namespace Ariadne {
  namespace Numeric {
    /*!
     * \brief A templated class representing an interval of real numbers.
     *
     * An interval of real numbers with endpoints of type \a R.
     * All operations on an interval must be guarenteed to return an interval contining the exact result.
     * If \a val of real numbers with endpoints of type \a R.
     * All operations on an interval must be guarenteed to return an interval contining the exact result.
     * If \a T supports exact evaluation of a function, then the exact evaluation must be used.
     * If \a T is dense in the reals, e.g. dyadic or rational, then any approximate operations may be given a maximum error of computation.
     *
     * Currently implemented as a wrapper around the boost::numeric::interval class template from the Boost C++ library.
     */
    template<typename R>
    class Interval : public boost::numeric::interval<R> {
     public:
      /*! \brief Default constructer constructs what?? */
      Interval()
        : boost::numeric::interval<R>(0,0) { }
      /*! \brief Construct a one-point interval. */
      Interval(const R& x)
        : boost::numeric::interval<R>(x) { }
      /*! \brief Construct from a boost interval. */
      Interval(const boost::numeric::interval<R>& ivl)
        : boost::numeric::interval<R>(ivl) { }
      /*! \brief Construct from lower and upper bounds. */
      Interval(const R& l, const R& u)
        : boost::numeric::interval<R>(l,u) { }
      /*! \brief Assignment operator. */
      const Interval<R>& operator=(const R& x) {
        *this=Interval<R>(x,x);
        return *this;
      }
      
      /*! \brief Equality operator. */
      bool operator==(const Interval<R>& ivl) const { 
        return this->lower()==ivl.lower() && this->upper()==ivl.upper(); }
      /*! \brief Inequality operator. */
      bool operator!=(const Interval<R>& ivl) const { 
        return !(*this==ivl); }

      /*! \brief Tests if the interval contains \a r. */
      bool contains(const R& r) const { return this->lower()<=r && r<=this->upper(); }
      /*! \brief Tests if the interior of the interval contains \a r. */
      bool interior_contains(const R& r) const { return this->lower()<r && r<this->upper(); }
      
      /*! \brief The midpoint of the interval, given by \f$(a+b)/2\f$. */
      R centre() const { return (this->lower()+this->upper())/2; }
      /*! \brief The radius of the interval, given by \f$(b-a)/2\f$. */
      R radius() const { return (this->upper()-this->lower())/2; }
      /*! \brief The length of the interval, given by \f$b-a\f$. */
      R length() const { return this->upper()-this->lower(); }
    };
    
    /*! \brief Tests disjointness. */
    template<typename R>
    inline
    bool
    disjoint(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.upper()<ivl2.lower() || ivl1.lower()>ivl2.upper());
    }

    /*! \brief Tests intersection of interiors. */
    template<typename R>
    inline
    bool
    interiors_intersect(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.upper()>ivl2.lower() && ivl1.lower()<ivl2.upper());
    }

    /*! \brief Tests if \a ivl1 is a subset of the interior of \a ivl2. */
    template<typename R>
    inline
    bool
    inner_subset(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.lower()>ivl2.lower() && ivl1.upper()<ivl2.upper());
    }

    /*! \brief Tests if \a ivl1 is a subset of \a ivl2. */
    template<typename R>
    inline
    bool
    subset(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.lower()>=ivl2.lower() && ivl1.upper()<=ivl2.upper());
    }

    /*! \brief Integer power. */
    template<typename R, typename N>
    inline
    Interval<R> 
    pow(const Interval<R>& x, const N& n) {
      Interval<R> result=R(1);
      for(N i=0; i!=n; ++i) {
        result*=x;
      }
      return result;
    }
  


    
    template<typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const Interval<R>& ivl)
    {
      const boost::numeric::interval<R>& boost_ivl(ivl);
      return os << boost_ivl;
    }
    
    template<typename R>
    inline
    std::istream&
    operator>>(std::istream& is, Interval<R>& ivl)
    {
      boost::numeric::interval<R>& boost_ivl(ivl);
      is >> boost_ivl;
      return is;
    }
    
  } // namespace Numeric
} // namespace Ariadne
  
#endif /* _ARIADNE_INTERVAL_H */
