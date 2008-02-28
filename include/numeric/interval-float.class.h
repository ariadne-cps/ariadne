/***************************************************************************
 *            interval.class.h
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
 
/*! \file interval.class.h
 *  \brief The class representing intervals of real number types. This header can be included to obtain declarations without needing to include definitions.
 */
 
#ifndef ARIADNE_INTERVAL_CLASS_H
#define ARIADNE_INTERVAL_CLASS_H

#include <iostream>
#include <stdexcept>

#include "base/tribool.h"
#include "numeric/expression.h"

namespace Ariadne {
  namespace Numeric {

  template<class R> class Interval;

    /*!\ingroup Numeric
     * \brief A templated class representing an interval of real numbers.
     * 
     * An interval of real numbers with endpoints of type \a R.
     * All operations on an interval must be guarenteed to return an interval contining the exact result.
     * If \a val of real numbers with endpoints of type \a R.
     * All operations on an interval must be guarenteed to return an interval contining the exact result.
     * If \a T supports exact evaluation of a function, then the exact evaluation must be used.
     * If \a T is dense in the reals, e.g. dyadic or rational, then any approximate operations may be given a max_imum error of computation.
     *
     * Currently implemented as a wrapper around the boost::numeric::interval class template from the Boost C++ library.
     */
    template<class T>
    class Interval< Float<T> >
      : public Value< Interval< Float<T> > >
    {
     private:
      typedef Float<T> R;
     public:
      R _lower; R _upper;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructer constructs empty interval. */
      Interval();
      /*! \brief Construct from a string literal. */
      Interval(const std::string& str);
      /*! \brief Construct a one-point interval. */
      Interval(const R& x);
      /*! \brief Construct from lower and upper bounds. */
      Interval(const R& l, const R& u);
      /*! \brief Copy constructor. */
      Interval(const Interval<R>& ivl);

      /*! \brief Assign from a number. */
      Interval<R>& operator=(const R& x);
      /*! \brief Copy assignment operator. */
      Interval<R>& operator=(const Interval<R>& ivl);

      /*! \brief Construct from an expression. */
      template<class E> Interval(const Expression<E>& e);
      /*! \brief Assign from an expression. */
      template<class E> Interval<R>& operator=(const Expression<E>& e);

      /*! \brief Construct a one-point interval. */
      template<class RX> Interval(const RX& x);
      /*! \brief Construct from lower and upper bounds. */
      template<class RL,class RU> Interval(const RL& l, const RU& u);
       /*! \brief Construct an interval with possibly different real type. */
      template<class RX> Interval(const Interval<RX>& ivl);
      /*! \brief Assign from a point of a possibly different type. */
      template<class RX> Interval<R>& operator=(const RX& x);
      /*! \brief Assign from an interval of possibly a different type. */
      template<class RX> Interval<R>& operator=(const Interval<RX>& ivl);
      //@}
      
      //@{
      //! \name Data access
      /*! \brief The lower bound. */
      const R& lower() const;
      /*! \brief The upper bound. */
      const R& upper() const;
      /*! \brief The midpoint of the interval, given by \f$(a+b)/2\f$. */
      R midpoint() const;
      /*! \brief The radius of the interval, given by \f$(b-a)/2\f$. */
      R radius() const;
      /*! \brief The width of the interval, given by \f$b-a\f$. */
      R width() const;
      //@}
      
      //@{
      //! \name Geometric operations
      
      /*! \brief Tests if the interval is empty. */
      bool empty() const;
      /*! \brief Tests if the interval consists of a single point. */
      bool singleton() const;
      /*! \brief Tests if the interval contains \a x. */
      template<class RX> bool encloses(const RX& x) const;
      /*! \brief Tests if the interval contains \a r. */
      template<class RX> bool refines(const Interval<RX>& ivl) const;
      //}
      
#ifdef DOXYGEN
      //@{
      //@{
      //! \name Approximation operations.
      /*! \brief The lower bound of the interval \a ivl. */
      friend R lower(const Interval<R>& ivl);
      /*! \brief The upper bound of the interval \a ivl. */
      friend R upper(const Interval<R>& ivl);
      /*! \brief The approximate midpoint of the interval \a ivl. */
      friend R midpoint(const Interval<R>& ivl);
      /*! \brief An upper bound for the error of the approximate value represented by the interval \a ivl. */
      friend R radius(const Interval<R>& ivl);
      /*! \brief An upper bound for the error of the approximate value represented by the interval \a ivl. */
      friend R width(const Interval<R>& ivl);
      /*! \brief Test if the interval \a ivl encloses the value \a x. */
      friend bool encloses(const Interval<R>& ivl, const R& x);
      /*! \brief Test if the interval \a ivl1 refines the approximation of the interval \a ivl2. */
      friend bool refines(const Interval<R>& ivl1, const Interval<R>& ivl2);
      //@}
      
      //! \name Arithmetic operations
      /*! \brief The interval of possible min_ima of \a x1 in \a ivl1 and \a x2 in \a ivl2. */
      friend Interval<R> min_(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The interval of possible max_ima of \a x1 in \a ivl1 and \a x2 in \a ivl2. */
      friend Interval<R> max_(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The interval of possible abs_olute values of \a x in \a ivl. */
      friend Interval<R> abs_(const Interval<R>& ivl);
      
      /*! \brief In-place addition of an interval. */
      friend Interval<R>& operator+=<>(Interval<R>&, const Interval<R>&);
      /*! \brief In-place addition of a number. */
      friend Interval<R>& operator+=<>(Interval<R>&, const R&);
      /*! \brief In-place subtraction of an interval. */
      friend Interval<R>& operator-=<>(Interval<R>&, const Interval<R>&);
      /*! \brief In-place subtraction of a number. */
      friend Interval<R>& operator-=<>(Interval<R>&, const R&);

      /*! \brief %Interval negation. */
      friend Interval<R> operator-(const Interval<R>& ivl);
      /*! \brief %Interval addition. */
      friend Interval<R> operator+(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief %Interval subtraction. */
      friend Interval<R> operator-(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief %Interval multiplication. */
      friend Interval<R> operator*(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief %Interval division. */
      friend Interval<R> operator/(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief %Integer power. */
      friend template<class N> Interval<R> pow_(const Interval<R>& x, const N& n);
      //@}
      
      //@{
      //! \name Geometric operations
      /*! \brief Tests equality. */
      friend bool equal(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests if \a ivl1 and \a ivl2 are disjoint. */
      friend bool disjoint(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests if the intersios of \a ivl1 and \a ivl2 overlap. */
      friend bool overlap(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests if \a ivl1 is a sub_set of \a ivl2. */
      friend bool subset(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests if \a ivl1 is a sub_set of the interior of \a ivl2. */
      friend bool inside(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests intersection of interiors. (%Deprecated) */
      friend bool interiors_intersect(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests if \a ivl1 is a sub_set of the interior of \a ivl2. (%Deprecated) */
      friend bool inner_subset(const Interval<R>& ivl1, const Interval<R>& ivl2);
      
      /*! \brief The intersection of \a ivl1 and \a ivl2. */
      friend Interval<R> intersection(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The closure of the intersection of the interiors of \a ivl1 and \a ivl2. (%Deprecated) */
      friend Interval<R> regular_intersection(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The smallest interval containing \a ivl1 and \a ivl2. */
      friend Interval<R> hull(const Interval<R>& ivl1, const Interval<R>& ivl2);
      //@}
      
      //@{
      //! \name Comparison operators.
      /*! \brief Equality operator. */
      friend template<class R1, class R2> tribool operator==(const Interval<R1>& ivl1, const Interval<R2>& ivl2); 
      /*! \brief Inequality operator. */
      friend template<class R1, class R2> tribool operator!=(const Interval<R1>& ivl1, const Interval<R2>& ivl2); 
      /*! \brief Less than operator. */
      friend template<class R1, class R2> tribool operator< (const Interval<R1>& ivl1, const Interval<R2>& ivl2);  
      /*! \brief Greater than operator. */
      friend template<class R1, class R2> tribool operator> (const Interval<R1>& ivl1, const Interval<R2>& ivl2);
      /*! \brief Less than or equal to operator. */
      friend template<class R1, class R2> tribool operator<=(const Interval<R1>& ivl1, const Interval<R2>& ivl2);
      /*! \brief Greater than or equal to operator. */
      friend template<class R1, class R2> tribool operator>=(const Interval<R1>& ivl1, const Interval<R2>& ivl2);

      /*! \brief Equality operator. */
      friend template<class R1, class R2> tribool operator==(const Interval<R1>& ivl1, const R2& ivl2); 
      /*! \brief Inequality operator. */
      friend template<class R1, class R2> tribool operator!=(const Interval<R1>& ivl1, const R2& ivl2); 
      /*! \brief Less than operator. */
      friend template<class R1, class R2> tribool operator< (const Interval<R1>& ivl1, const R2& ivl2);  
      /*! \brief Greater than operator. */
      friend template<class R1, class R2> tribool operator> (const Interval<R1>& ivl1, const R2& ivl2);
      /*! \brief Less than or equal to operator. */
      friend template<class R1, class R2> tribool operator<=(const Interval<R1>& ivl1, const R2& ivl2);
      /*! \brief Greater than or equal to operator. */
      friend template<class R1, class R2> tribool operator>=(const Interval<R1>& ivl1, const R2& ivl2);

      /*! \brief Equality operator. */
      friend template<class R1, class R2> tribool operator==(const R1& ivl1, const Interval<R2>& ivl2); 
      /*! \brief Inequality operator. */
      friend template<class R1, class R2> tribool operator!=(const R1& ivl1, const Interval<R2>& ivl2); 
      /*! \brief Less than operator. */
      friend template<class R1, class R2> tribool operator< (const R1& ivl1, const Interval<R2>& ivl2);  
      /*! \brief Greater than operator. */
      friend template<class R1, class R2> tribool operator> (const R1& ivl1, const Interval<R2>& ivl2);
      /*! \brief Less than or equal to operator. */
      friend template<class R1, class R2> tribool operator<=(const R1& ivl1, const Interval<R2>& ivl2);
      /*! \brief Greater than or equal to operator. */
      friend template<class R1, class R2> tribool operator>=(const R1& ivl1, const Interval<R2>& ivl2);
      //@}

      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const Interval<R>& ivl);
      /*! \brief Stream extraction operator. */
      friend std::istream& operator>>(std::istream& is, Interval<R>& ivl);
      //@}
#endif
    };

    template<class X> std::string name();

    


  } 
}


#endif /* ARIADNE_INTERVAL_CLASS_H */
