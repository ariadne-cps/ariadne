/***************************************************************************
 *            approximation.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file numerical_type.h
 *  \brief Type definitions and conversion operators for fundamental Ariadne types.
 */

#ifndef _ARIADNE_APPROXIMATION_H
#define _ARIADNE_APPROXIMATION_H

#include <gmpxx.h>
#include <string>
#include <assert.h>

#include "../base/numerical_type.h"

namespace Ariadne {

  /*! \brief Approximate \a x by an element of \p Res with accuracy \a e. */
  template<typename Res, typename Arg, typename Err> Res approximate(const Arg& x, const Err& e);
  
  template<> inline double approximate(const double& x, const double& e) {
    return x;
  }
  
  template<> inline Rational approximate(const Rational& q, const Rational& e) {
    return q;
  }
  
  template<> inline Dyadic approximate(const Rational& q, const Rational& e) {
    Dyadic x=Dyadic(q);
    assert(abs(Rational(x)-q)<e);
    return x;
  }
    
  template<> inline Dyadic approximate(const Rational& q, const Dyadic& e) {
    return approximate<Dyadic>(q,Rational(e));
  }
    
  template<> inline MPFloat approximate(const Rational& q, const Rational& e) {
    MPFloat x=MPFloat(q);
    assert(abs(Rational(x)-q)<e);
    return x;
  }
    
  template<> inline MPFloat approximate(const Rational& q, const MPFloat& e) {
    return approximate<MPFloat>(q,Rational(e));
  }
    
}

#endif /* _ARIADNE_APPROXIMATION_TYPE */
