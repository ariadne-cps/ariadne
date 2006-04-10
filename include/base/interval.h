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
  namespace Base {
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
     * COMMENT FOR ALBERTO: How should we specify error bounds for computations on types which support arbitrary-precision computing?
     * I would suggest using a global(ish) "precision" variable which either may be "locked" by a computation, or stored on a stack.
     * I think functions should keep a natural syntax, so we can declare <code>rational cos(rational)</code>,
     * and set error bounds for the computation elsewhere.
     * Of course, this partly depends on the interval / rational arithmetic library we use.
     *
     * Currently implemented as a wrapper around the boost::numeric::interval class template from the Boost C++ library.
     */
    template<typename R>
    class Interval : public boost::numeric::interval<R> {
     public:
      Interval()
        : boost::numeric::interval<R>() { }
      Interval(const boost::numeric::interval<R>& ivl)
        : boost::numeric::interval<R>(ivl) { }
      Interval(const R& l, const R& u)
        : boost::numeric::interval<R>(l,u) { }
      
      bool operator==(const Interval<R>& ivl) { 
        return this->lower()==ivl.lower() && this->upper()==ivl.upper(); }
      bool operator!=(const Interval<R>& ivl) { 
        return !(*this==ivl); }

      const Interval<R>& operator=(const R &r) {
        *this=Interval<R>(r,r);
        return *this;
      }

      bool contains(const R& r) const { return this->lower()<=r && r<=this->upper(); }
      
      R centre() const { return (this->lower()+this->upper())/2; }
      R radius() const { return (this->upper()-this->lower())/2; }
      R length() const { return this->upper()-this->lower(); }
    };
    
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
    
  } // namespace Base
} // namespace Ariadne
  
#endif /* _ARIADNE_INTERVAL_H */
