/***************************************************************************
 *            interval.h
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
 
/*! \file interval.h
 *  \brief Intervals representing sets of real numbers.
 */
 
#ifndef ARIADNE_GEOMETRY_INTERVAL_H
#define ARIADNE_GEOMETRY_INTERVAL_H

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cassert>

#include "../base/tribool.h"
#include "../base/exceptions.h"


namespace Ariadne {
  namespace Geometry {
  
 
    /*!\ingroup Geometry
     * \brief A templated class representing an interval set of real numbers. 
     * 
     * The interval is specified by the lower and upper endpoints, and must be binary compatible with the Numeric::FuzzyFloat class.
     */
    template<class X>
    class Interval
    {
     private:
      X _lower; X _upper;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructer constructs empty interval. */
      Interval();
      /*! \brief Construct a singleton set based on the point \a x. */
      Interval(const X& x);
      /*! \brief Construct from lower and upper bounds. */
      Interval(const X& l, const X& u);
      /*! \brief Copy constructor. */
      Interval(const Interval<X>& ivl);
      /*! \brief Assign from a number. */
      Interval<X>& operator=(const X& x);
      /*! \brief Copy assignment operator. */
      Interval<X>& operator=(const Interval<X>& ivl);
 
      /*! \brief Construct a one-point interval. */
      template<class XR> Interval(const XR& x);
      /*! \brief Construct from lower and upper bounds. */
      template<class XL,class XU> Interval(const XL& l, const XU& u);
       /*! \brief Construct an interval with possibly different real type. */
      template<class XR> Interval(const Interval<XR>& ivl);
      /*! \brief Assign from a value of a possibly different type. */
      template<class XR> Interval<X>& operator=(const XR& x);
      /*! \brief Assign from an interval of possibly a different type. */
      template<class XR> Interval<X>& operator=(const Interval<XR>& ivl);
      //@}
      
      //@{
      //! \name Comparison operators
      /*! \brief Equality operator tests exact equality of representation. */
      bool operator==(const Interval<X>& ivl);
      /*! \brief Inequality operator compares equality of representation. */
      bool operator!=(const Interval<X>& ivl);
      //@}
      
      //@{
      //! \name Data access
      /*! \brief The lower bound. */
      const X& lower() const;
      /*! \brief The lower bound. */
      const X& lower_bound() const;
      /*! \brief The upper bound. */
      const X& upper() const;
      /*! \brief The upper bound. */
      const X& upper_bound() const;
      //@}
      
      //@{
      //! \name Input/output operations
      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);
      //@}

#ifdef DOXYGEN
     
      //@{
      //! \name Geometric operations
      /*! \brief The approximate centre of the interval \a ivl. */
      friend X centre<>(const Interval<X>& ivl);
      /*! \brief An upper bound for the radius of \a ivl. */
      friend X radius<>(const Interval<X>& ivl);
      /*! \brief Tests if \a ivl is empty. */
      friend tribool empty<>(const Interval<X>& ivl);
      /*! \brief Tests if \a ivl is bounded. */
      friend tribool bounded<>(const Interval<X>& ivl);
      /*! \brief Tests if \a ivl1 contains the point \a x. */
      friend tribool contains<>(const Interval<X>& ivl, const X& x);
      /*! \brief Tests if \a ivl1 and \a ivl2 are disjoint. */
      friend tribool disjoint<>(const Interval<X>& ivl1, const Interval<X>& ivl2);
      /*! \brief Tests if \a ivl1 and \a ivl2 intersect. */
      friend tribool intersect<>(const Interval<X>& ivl1, const Interval<X>& ivl2);
      /*! \brief Tests if \a ivl1 is a subset of \a ivl2. */
      friend tribool subset<>(const Interval<X>& ivl1, const Interval<X>& ivl2);
      /*! \brief Tests if \a ivl1 is a superset of \a ivl2. */
      friend tribool superset<>(const Interval<X>& ivl1, const Interval<X>& ivl2);
      
      /*! \brief The intersection of \a ivl1 and \a ivl2. */
      friend Interval<X> closed_intersection<>(const Interval<X>& ivl1, const Interval<X>& ivl2);
      /*! \brief The closure of the intersection of the interiors of \a ivl1 and \a ivl2. (%Deprecated) */
      friend Interval<X> open_intersection<>(const Interval<X>& ivl1, const Interval<X>& ivl2);
      /*! \brief The smallest interval containing \a ivl1 and \a ivl2. */
      friend Interval<X> interval_hull<>(const Interval<X>& ivl1, const Interval<X>& ivl2);
      //@}
      
      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const Interval<X>& ivl);
      /*! \brief Stream extraction operator. */
      friend std::istream& operator>>(std::istream& is, Interval<X>& ivl);
      //@}
#endif
    };
    

    template<class X> X lower(const Interval<X>& x);
    template<class X> X upper(const Interval<X>& x);
    template<class X> X centre(const Interval<X>& x);
    template<class X> X radius(const Interval<X>& x);

    template<class X> tribool empty(const Interval<X>& x);
    template<class X> tribool bounded(const Interval<X>& x);

    template<class X1, class X2> tribool contains(const Interval<X1>& x1, const X2& x2);
    template<class X1, class X2> tribool disjoint(const Interval<X1>& x1, const Interval<X2>& x2);
    template<class X1, class X2> tribool intersect(const Interval<X1>& x1, const Interval<X2>& x2);
    template<class X1, class X2> tribool subset(const Interval<X1>& x1, const Interval<X2>& x2);
    template<class X1, class X2> tribool superset(const Interval<X1>& x1, const Interval<X2>& x2);
    template<class X> Interval<X> closed_intersection(const Interval<X>& x1, const Interval<X>& x2);
    template<class X> Interval<X> open_intersection(const Interval<X>& x1, const Interval<X>& x2);
    template<class X> Interval<X> interval_hull(const Interval<X>& x1, const Interval<X>& x2);
    
    template<class X> std::ostream& operator<<(std::ostream& os, const Interval<X>& x);
    template<class X> std::istream& operator>>(std::istream& is, Interval<X>& x);
    
  } // namespace Geometry
} // namespace Ariadne

#include "interval.inline.h"

#endif /* ARIADNE_GEOMETRY_INTERVAL_H */
