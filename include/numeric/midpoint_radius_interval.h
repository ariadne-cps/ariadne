/***************************************************************************
 *            midpoint_radius_interval.h
 *
 *  Copyright 2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file midpoint_radius_interval.h
 *  \brief Intervals of real number types using midpoint-radius representations.
 */
 
#ifndef ARIADNE_MIDPOINT_RADIUS_INTERVAL_H
#define ARIADNE_MIDPOINT_RADIUS_INTERVAL_H

#include <iostream>
#include <stdexcept>

#include "../base/tribool.h"
#include "../base/exceptions.h"

#include "../numeric/exceptions.h"
#include "../numeric/traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"

namespace Ariadne {
  namespace Numeric {

    /*!\ingroup Numeric
     * \brief A templated class representing an interval of real numbers using the midpoint-radius representation.
     * 
     * An interval of real numbers with midpoint and radius of type \a R.
     * All operations on an interval must be guarenteed to return an interval contining the exact result.
     * If \a T supports exact evaluation of a function, then the exact evaluation must be used.
     * If \a T is dense in the reals, e.g. dyadic or rational, then any approximate operations may be given a max_imum error of computation.
     */
    template<class R>
    class MidpointRadiusInterval
    {
     private:
      R _midpoint; R _radius;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructer constructs empty interval. */
      MidpointRadiusInterval() : _midpoint(0), _radius(0) { }
      /*! \brief Default constructer constructs empty interval. */
      MidpointRadiusInterval(const MidpointRadiusInterval<R>& ivl) : _midpoint(ivl._midpoint), _radius(ivl._radius) { }
      /*! \brief Construct from a single point. */
      MidpointRadiusInterval(const R& x) : _midpoint(x), _radius(0) { }
      /*! \brief Construct from a midpoint and a radius. */
      MidpointRadiusInterval(const R& x, const R& r) : _midpoint(x), _radius(r) { }
      /*! \brief Default constructer constructs empty interval. */
      MidpointRadiusInterval(const Interval<R>& ivl) : _midpoint(ivl.midpoint()), _radius(ivl.radius()) { }
      /*! \brief Construct from a single point. */
      MidpointRadiusInterval& operator=(const R& x) { this->_midpoint=x; this->_radius=0; return *this; }
      /*! \brief Default constructer constructs empty interval. */
      MidpointRadiusInterval& operator=(const Interval<R>& ivl) { this->_midpoint=ivl.midpoint(); this->_radius=ivl.radius(); return *this; }
      //@}
      
      //@{
      //! \name Data access
      /*! \brief The lower bound. */
      const R& midpoint() const { return this->_midpoint; }
      /*! \brief The upper bound. */
      const R& radius() const { return this->_radius; }
      /*! \brief The lower bound. */
      R lower() const { return sub_down(this->_midpoint,this->_radius); }
      /*! \brief The upper bound. */
      R upper() const { return add_up(this->_midpoint,this->_radius); }
      //@}
      
#ifdef DOXYGEN
      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const Interval<R>& ivl);
      /*! \brief Stream extraction operator. */
      friend std::istream& operator>>(std::istream& is, Interval<R>& ivl);
      //@}
#endif
    };
 
    
    template<class R>
    MidpointRadiusInterval<R> dot(int n, const R* ptrX, int incX, const R* ptrY, int incY) {
      Interval<R> r;
      for(int i=0; i!=n; ++i) {
        r+=ptrX[i*n]+ptrY[i*n];
      }
      return MidpointRadiusInterval<R>(r);
    }

   
  }
}
#endif /* ARIADNE_MIDPOINT_RADIUS_INTERVAL_H */
