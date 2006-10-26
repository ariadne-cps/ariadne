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
 *  \brief Simple arithmetic functions for integer, dyadic and rational types and intervals.
 */
 
#ifndef _ARIADNE_ARITHMETIC_H
#define _ARIADNE_ARITHMETIC_H

#include <algorithm>
#include <cassert>

namespace Ariadne {
  namespace Numeric {

    template<class R> class Interval;

    //! \name Exact arithmetical operations. 
    //@{
    //! \ingroup Numeric
    /*! \brief Minimum. */
    template<class R> inline R min(const R& x1, const R& x2);
    template<class R> inline Interval<R> min(const Interval<R>& x1, 
                                                const Interval<R>& x2);

    /*! \brief Maximum. */
    template<class R> inline R max(const R& x1, const R& x2);;
    template<class R> inline Interval<R> max(const Interval<R>& x1, 
                                                const Interval<R>& x2);

    /*! \brief The median (average) of two values. */
    template<class R> inline R med(const R& x, const R& y);

    /*! \brief Absolute value. */
    template<class R> inline R abs(const R& x);  
    template<class R> inline Interval<R> abs(const Interval<R>& x1);

    /*! \brief Unary negation. */
    template<class R> inline R neg(const R& x);
    
    /*! \brief Addition. */
    template<class R> inline R add(const R& x1,const R& x2);
    
    /*! \brief Subtraction. */
    template<class R> inline R sub(const R& x1,const R& x2);
    
    /*! \brief Multiplication. */
    template<class R> inline R mul(const R& x1,const R& x2);
    
    /*! \brief Division. */
    template<class R> inline R div(const R& x1, const R& x2);
    
    /*! \brief The power of a real number type by an integer. */
    template<class R, class N> inline R pow(const R& x, const N& n);
    
    /*! \brief The integer part of \a x, rounding towards 0. */
    template<class N, class R> inline N floor(const R& x);

    /*! \brief The integer ceiling of \a x, rounding away from 0. */
    template<class N, class R> inline N ceil(const R& x);
    
    /*! \brief An integer \a n such that \f$n\leq x/y < n+1\f$. */
    template<class N, class R> inline N quot(const R& x, const R& y);
    //@}

  }
}



#endif /* _ARITHMETIC_H */
