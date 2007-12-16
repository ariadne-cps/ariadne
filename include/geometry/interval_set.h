/***************************************************************************
 *            geometry/interval_set.h
 *
 *  Copyright 2007  Pieter Collins
 *
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
 
/*! \file geometry/interval_set.h
 *  \brief Intervals representing sets of real numbers.
 */
 
#ifndef ARIADNE_GEOMETRY_INTERVAL_SET_H
#define ARIADNE_GEOMETRY_INTERVAL_SET_H

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cassert>

#include "base/tribool.h"
#include "base/exceptions.h"



namespace Ariadne {
  namespace Geometry {
  
 
    /*!\ingroup Geometry
     * \brief A templated class representing an interval set of real numbers. 
     * 
     * The interval is specified by the lower and upper endpoints, and must be binary compatible with the Numeric::FuzzyFloat class.
     */
    template<class X>
    class IntervalSet
    {
     private:
      X _lower; X _upper;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructer constructs empty interval. */
      IntervalSet();
      /*! \brief Construct a singleton set based on the point \a x. */
      IntervalSet(const X& x);
      /*! \brief Construct from lower and upper bounds. */
      IntervalSet(const X& l, const X& u);
      /*! \brief Copy constructor. */
      IntervalSet(const IntervalSet<X>& ivl);
      /*! \brief Assign from a number. */
      IntervalSet<X>& operator=(const X& x);
      /*! \brief Copy assignment operator. */
      IntervalSet<X>& operator=(const IntervalSet<X>& ivl);
 
      /*! \brief Construct a one-point interval. */
      template<class XR> IntervalSet(const XR& x);
      /*! \brief Construct from lower and upper bounds. */
      template<class XL,class XU> IntervalSet(const XL& l, const XU& u);
       /*! \brief Construct an interval with possibly different real type. */
      template<class XR> IntervalSet(const IntervalSet<XR>& ivl);
      /*! \brief Assign from a value of a possibly different type. */
      template<class XR> IntervalSet<X>& operator=(const XR& x);
      /*! \brief Assign from an interval of possibly a different type. */
      template<class XR> IntervalSet<X>& operator=(const IntervalSet<XR>& ivl);
      //@}
      
      //@{
      //! \name Comparison operators
      /*! \brief Equality operator tests exact equality of representation. */
      bool operator==(const IntervalSet<X>& ivl);
      /*! \brief Inequality operator compares equality of representation. */
      bool operator!=(const IntervalSet<X>& ivl);
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
      friend X centre<>(const IntervalSet<X>& ivl);
      /*! \brief An upper bound for the radius of \a ivl. */
      friend X radius<>(const IntervalSet<X>& ivl);
      /*! \brief Tests if \a ivl is empty. */
      friend tribool empty<>(const IntervalSet<X>& ivl);
      /*! \brief Tests if \a ivl is bounded. */
      friend tribool bounded<>(const IntervalSet<X>& ivl);
      /*! \brief Tests if \a ivl1 contains the point \a x. */
      friend tribool contains<>(const IntervalSet<X>& ivl, const X& x);
      /*! \brief Tests if \a ivl1 and \a ivl2 are disjoint. */
      friend tribool disjoint<>(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2);
      /*! \brief Tests if \a ivl1 and \a ivl2 intersect. */
      friend tribool intersect<>(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2);
      /*! \brief Tests if \a ivl1 is a subset of \a ivl2. */
      friend tribool subset<>(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2);
      /*! \brief Tests if \a ivl1 is a superset of \a ivl2. */
      friend tribool superset<>(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2);
      
      /*! \brief The intersection of \a ivl1 and \a ivl2. */
      friend IntervalSet<X> closed_intersection<>(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2);
      /*! \brief The closure of the intersection of the interiors of \a ivl1 and \a ivl2. (%Deprecated) */
      friend IntervalSet<X> open_intersection<>(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2);
      /*! \brief The smallest interval containing \a ivl1 and \a ivl2. */
      friend IntervalSet<X> interval_hull<>(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2);
      //@}
      
      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const IntervalSet<X>& ivl);
      /*! \brief Stream extraction operator. */
      friend std::istream& operator>>(std::istream& is, IntervalSet<X>& ivl);
      //@}
#endif
    };
    

    template<class X> X lower(const IntervalSet<X>& x);
    template<class X> X upper(const IntervalSet<X>& x);
    template<class X> X centre(const IntervalSet<X>& x);
    template<class X> X radius(const IntervalSet<X>& x);

    template<class X> tribool empty(const IntervalSet<X>& x);
    template<class X> tribool bounded(const IntervalSet<X>& x);

    template<class X1, class X2> tribool contains(const IntervalSet<X1>& x1, const X2& x2);
    template<class X1, class X2> tribool disjoint(const IntervalSet<X1>& x1, const IntervalSet<X2>& x2);
    template<class X1, class X2> tribool intersect(const IntervalSet<X1>& x1, const IntervalSet<X2>& x2);
    template<class X1, class X2> tribool subset(const IntervalSet<X1>& x1, const IntervalSet<X2>& x2);
    template<class X1, class X2> tribool superset(const IntervalSet<X1>& x1, const IntervalSet<X2>& x2);
    template<class X> IntervalSet<X> closed_intersection(const IntervalSet<X>& x1, const IntervalSet<X>& x2);
    template<class X> IntervalSet<X> open_intersection(const IntervalSet<X>& x1, const IntervalSet<X>& x2);
    template<class X> IntervalSet<X> interval_hull(const IntervalSet<X>& x1, const IntervalSet<X>& x2);
    
    template<class X> std::ostream& operator<<(std::ostream& os, const IntervalSet<X>& x);
    template<class X> std::istream& operator>>(std::istream& is, IntervalSet<X>& x);
    
  } // namespace Geometry
} // namespace Ariadne

#include "interval_set.inline.h"

#endif /* ARIADNE_GEOMETRY_INTERVAL_SET_H */
