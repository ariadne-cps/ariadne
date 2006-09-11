/***************************************************************************
 *            arithmetic.h
 *
 *  Wed 18 May 2005
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
 
/*! \file arithmetic.h
 *  \brief Simple arithmetic functions for integer, double, dyadic and rational types and intervals.
 */
 
#ifndef _ARIADNE_ARITHMETIC_H
#define _ARIADNE_ARITHMETIC_H

#include <algorithm>
#include <cassert>

namespace Ariadne {
  namespace Numeric {

    template<typename R> class Interval;

    //! \name Exact arithmetical operations. 
    //@{
    //! \ingroup Numeric
    /*! \brief Minimum. */
    template<typename R, typename A1, typename A2>
    inline
    R 
    min(const A1& x1, const A1& x2) {
      //std::cerr << "min<" << name<R>() << ">" << std::endl;
      return x1<x2 ? R(x1) : R(x2);
    }
  
    template<typename R>
    inline
    R 
    min(const R& x1, const R& x2) {
      //std::cerr << "min<" << name<R>() << ">" << std::endl;
      return x1<x2 ? x1 : x2;
    }
  
    template<typename R>
    inline
    Interval<R> 
    min(const Interval<R>& x1, const Interval<R>& x2);

  
    /*! \brief Maximum. */
    template<typename R, typename A1, typename A2>
    inline
    R 
    max(const A1& x1, const A1& x2) {
      //std::cerr << "min<" << name<R>() << ">" << std::endl;
      return x1>x2 ? R(x1) : R(x2);
    }
  
    template<typename R>
    inline
    R 
    max(const R& x1, const R& x2) {
      //std::cerr << "max<" << name<R>() << ">" << std::endl;
      return x1>x2 ? x1 : x2;
    }
  
    template<typename R>
    inline
    Interval<R> 
    max(const Interval<R>& x1, const Interval<R>& x2);

 
    /*! \brief The median (average) of two values. */
    template<typename R>
    inline
    R 
    med(const R& x, const R& y) {
      return (x+y)/2;
    }

    /*! \brief Absolute value. */
    template<typename R>
    inline
    R 
    abs(const R& x) {
      return x>=R(0) ? x : R(-x);
    }
  
    template<typename R>
    inline
    Interval<R> 
    abs(const Interval<R>& x1, const Interval<R>& x2);


    /*! \brief Unary negation. */
    template<typename R> 
    inline
    R neg(const R& x) {
      return -x;
    }

    
    /*! \brief Addition. */
    template<typename R> 
    inline
    R add(const R& x1,const R& x2) {
      return x1+x2;
    }

    
    /*! \brief Subtraction. */
    template<typename R> 
    inline
    R sub(const R& x1,const R& x2) {
      return x1-x2;
    }

    
    /*! \brief Multiplication. */
    template<typename R> 
    inline
    R mul(const R& x1,const R& x2) {
      return x1*x2;
    }

    
    /*! \brief Division. */
    template<typename R> 
    inline
    R
    div(const R& x1, const R& x2)
    {
      return x1/x2;
    }

    /*! \brief The power of a real number type by an integer. */
    template<typename R, typename N> R pow(const R& x, const N& n) {
      R r=1; R p=x; N e=1;
      while(e<n) { if(e&n) { r*=p; } p*=p; e*=2; }
      return r; 
    }
    
    /*! \brief The integer part of \a x, rounding towards 0. */
    template<typename N, typename R> N floor(const R& x);

    /*! \brief The integer ceiling of \a x, rounding away from 0. */
    template<typename N, typename R> N ceil(const R& x);
    
    /*! \brief An integer \a n such that \f$n\leq x/y < n+1\f$. */
    template<typename N, typename R> N quot(const R& x, const R& y);

    /*! \brief An integer \a n such that \f$n\leq x/y < n+1\f$. (Deprecated) 
     *  \deprecated This function has been replaced by quot . 
     */
    template<typename R> inline int quotient(const R& x, const R& y) { return quot<int>(x,y); }
    

    //@}

  }
}



#endif /* _ARITHMETIC_H */
