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
 *  \brief Intervals of real number types (currently implemented using Boost).
 */
 
#ifndef ARIADNE_INTERVAL_H
#define ARIADNE_INTERVAL_H

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cassert>

#include "../base/tribool.h"
#include "../base/exceptions.h"

#include "../numeric/exceptions.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"

#include "../numeric/interval.class.h"

namespace Ariadne {
  namespace Numeric {
  
    /*!\ingroup Numeric
     * \brief A reference to an interval. 
     */
    template<class R>
    class IntervalReference {
     public:
      IntervalReference(R& l, R& u);
      IntervalReference(Interval<R>& ivl);
      void operator=(const Interval<R>& ivl);
      operator Interval<R> () const;
      R lower() const;
      R upper() const;
      R centre() const;
     private:
      R* _lower; R* _upper;
    };
    
    
  

    template<class R> R lower(const Interval<R>& x);
    template<class R> R upper(const Interval<R>& x);
    template<class R> R midpoint(const Interval<R>& x);
    template<class R> R radius(const Interval<R>& x);
    template<class R> R width(const Interval<R>& x);

    template<class R, class RX> bool encloses(const Interval<R>& ivl, const RX& x);
    template<class R1, class R2> bool refines(const Interval<R1>& ivl1, const Interval<R2>& ivl2);

    template<class R1, class R2> tribool operator==(const Interval<R1>& ivl1, const Interval<R2>& ivl2);
    template<class R1, class R2> tribool operator!=(const Interval<R1>& ivl1, const Interval<R2>& ivl2);
    template<class R1, class R2> tribool operator<(const Interval<R1>& ivl1, const Interval<R2>& ivl2);
    template<class R1, class R2> tribool operator>(const Interval<R1>& ivl1, const Interval<R2>& ivl2);
    template<class R1, class R2> tribool operator<=(const Interval<R1>& ivl1, const Interval<R2>& ivl2);
    template<class R1, class R2> tribool operator>=(const Interval<R1>& ivl1, const Interval<R2>& ivl2);

    template<class R1, class R2> tribool operator==(const Interval<R1>& ivl, const R2& x);
    template<class R1, class R2> tribool operator!=(const Interval<R1>& ivl, const R2& x);
    template<class R1, class R2> tribool operator<(const Interval<R1>& ivl, const R2& x);
    template<class R1, class R2> tribool operator>(const Interval<R1>& ivl, const R2& x);
    template<class R1, class R2> tribool operator<=(const Interval<R1>& ivl, const R2& x);
    template<class R1, class R2> tribool operator>=(const Interval<R1>& ivl, const R2& x);

    template<class R1, class R2> tribool operator==(const R1& x, const Interval<R2>& ivl);
    template<class R1, class R2> tribool operator!=(const R1& x, const Interval<R2>& ivl);
    template<class R1, class R2> tribool operator<(const R1& x, const Interval<R2>& ivl);
    template<class R1, class R2> tribool operator>(const R1& x, const Interval<R2>& ivl);
    template<class R1, class R2> tribool operator<=(const R1& x, const Interval<R2>& ivl);
    template<class R1, class R2> tribool operator>=(const R1& x, const Interval<R2>& ivl);

    template<class R> Interval<R> operator+(const Interval<R>& x);
    template<class R> Interval<R> operator-(const Interval<R>& x);
    
    template<class R> Interval<R> operator+(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator+(const Interval<R>& x1, const R& x2);
    template<class R> Interval<R> operator+(const R& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator+=(Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator+=(Interval<R>& x1, const R& x2);

    template<class R> Interval<R> operator-(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator-(const Interval<R>& x1, const R& x2);
    template<class R> Interval<R> operator-(const R& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator-=(Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator-=(Interval<R>& x1, const R& x2);

    template<class R> Interval<R> operator*(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator*(const Interval<R>& x1, const R& x2);
    template<class R> Interval<R> operator*(const R& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator*=(Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator*=(Interval<R>& x1, const R& x2);

    template<class R> Interval<R> operator/(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator/(const Interval<R>& x1, const R& x2);
    template<class R> Interval<R> operator/(const R& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator/=(Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator/=(Interval<R>& x1, const R& x2);

    template<class R> Interval<R> operator+(const Interval<R>& x1, const int& x2);
    template<class R> Interval<R> operator-(const Interval<R>& x1, const int& x2);
    template<class R> Interval<R> operator*(const Interval<R>& x1, const int& x2);
    template<class R> Interval<R> operator/(const Interval<R>& x1, const int& x2);

    template<class R> Interval<R> operator+(const int& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator-(const int& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator*(const int& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator/(const int& x1, const Interval<R>& x2);
    
    template<class R> Interval<R> operator+(const Interval<R>& x1, const double& x2);
    template<class R> Interval<R> operator-(const Interval<R>& x1, const double& x2);
    template<class R> Interval<R> operator*(const Interval<R>& x1, const double& x2);
    template<class R> Interval<R> operator/(const Interval<R>& x1, const double& x2);

    template<class R> Interval<R> operator+(const double& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator-(const double& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator*(const double& x1, const Interval<R>& x2);
    template<class R> Interval<R> operator/(const double& x1, const Interval<R>& x2);

    template<class R> Interval<R> min(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> max(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> abs(const Interval<R>& x);
    template<class R,class N> Interval<R> pow(const Interval<R>& x, const N& n);
  
    template<class R> bool equal(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> bool disjoint(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> bool overlap(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> bool subset(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> bool inside(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> intersection(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> hull(const Interval<R>& x1, const Interval<R>& x2);
    
    template<class R> Interval<R> sqrt(const Interval<R>& x);
    template<class R> Interval<R> hypot(const Interval<R>& x1, const Interval<R>& x2);
    template<class R> Interval<R> exp(const Interval<R>& x);
    template<class R> Interval<R> log(const Interval<R>& x);

    template<typename R> Interval<R> pi();
    template<typename R> Interval<R> sin(const Interval<R>& ivl);
    template<typename R> Interval<R> cos(const Interval<R>& ivl);
    template<typename R> Interval<R> tan(const Interval<R>& ivl);
    template<typename R> Interval<R> asin(const Interval<R>& ivl);
    template<typename R> Interval<R> acos(const Interval<R>& ivl);
    template<typename R> Interval<R> atan(const Interval<R>& ivl);
    
    template<class R> std::ostream& operator<<(std::ostream& os, const Interval<R>& x);
    template<class R> std::istream& operator>>(std::istream& is, Interval<R>& x);
    
  } // namespace Numeric
} // namespace Ariadne
  

namespace TBLAS {
  template<class real> int iamax (const int N, const real *X, const int incX);
  template<class real> int iamax (const int N, const Ariadne::Numeric::Interval<real> *X, const int incX);
}

#include "interval.inline.h"

#endif /* ARIADNE_INTERVAL_H */
