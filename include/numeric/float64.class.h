/***************************************************************************
 *            numeric/float64.h
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
 
/*! \file numeric/float64.h
 *  \brief Type definitions and conversion operators for 64-bit fixed precision floating point numbers.
 */

#ifndef ARIADNE_NUMERIC_FLOAT64_CLASS_H
#define ARIADNE_NUMERIC_FLOAT64_CLASS_H

#include <iosfwd>

#include <boost/numeric/interval/rounded_arith.hpp>
#include <boost/numeric/interval/rounded_transc.hpp>

// Try uncommenting this code if experiencing difficulties with rounded and interval transcendental functions.
/*
#ifdef __x86_64__ 
#ifndef __i386__
#define __i386__
#warning "Defining __i386__ preprocessor symbol in order to use x86 hardware rounding on an __x86_64__ machine. Please contact the developers if this causes a compilation problem."
#include <boost/numeric/interval/hw_rounding.hpp>
#undef __i386__
#endif
#endif
*/

#include <boost/numeric/interval/hw_rounding.hpp>

#include "numeric/expression.h"

namespace Ariadne {
  namespace Numeric {

    class Integer;
    class Rational;
    template<class T> class Float;
    template<class R> class Interval;

    typedef Float<double> Float64;
    typedef Interval<Float64> Interval64;


    /*!\ingroup Numeric
     * \brief A 64-bit fixed-precision floating point number.
     *  
     * Standard operations are not exact, but must support interval arithmetic.
     *
     * Currently implemented by the built-in type double.
     */
    template<>
    class Float<double> 
      : public Value< Float<double> >
    {
     public:
      typedef boost::numeric::interval_lib::rounded_arith_std<double> rounded_arith;
      typedef boost::numeric::interval_lib::rounded_transc_std<double> rounding;
     public:
      ~Float();
      Float();
      Float(const int& n);
      Float(const unsigned int& n);
      Float(const double& x);
      Float(const Integer& n);
      Float(const Float64& x);
      
      explicit Float(const std::string& x);
      
      Float64& operator=(const int& n);
      Float64& operator=(const unsigned int& n);
      Float64& operator=(const double& x);
      Float64& operator=(const Integer& n);
      Float64& operator=(const Float64& x);
      
      template<class E> Float(const Expression<E>& e);
      template<class E> Float64& operator=(const Expression<E>& e);

      template<class X, class Rnd> Float(const X& x, Rnd rnd);
      template<class E, class Rnd> Float(const Expression<E>& x, Rnd rnd);

      double get_d() const;
     public:
      double _value;
    };

    std::ostream& operator<<(std::ostream& os, const Float64& x);
    std::istream& operator>>(std::istream& is, Float64& x);

  }
}

#endif /* ARIADNE_NUMERIC_FLOAT64_CLASS_H */
