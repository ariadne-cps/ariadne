/***************************************************************************
 *            integer.h
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
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
 
/*! \file integer.h
 *  \brief Multiple-precision integer type and interger functions.
 */

#ifndef _ARIADNE_INTEGER_H
#define _ARIADNE_INTEGER_H

#include <gmpxx.h>

#include "../declarations.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"

namespace Ariadne {
  namespace Numeric {

#ifdef DOXYGEN
    /*!\ingroup Numeric
     * \brief An integer of arbitrary size.
     *  
     * An element of the ring of integers.
     * Must allow denotation of any integer, including arbitrarily large values.
     * Integer quotient and remainder must be supported.
     *
     * Currently implemented using mpz_class from the GNU Multiple Precision Library.
     */
    class Integer { };
#else
    typedef mpz_class Integer;
#endif
  
    template<> class numerical_traits<Integer> {
     public:
      typedef ring_tag algebraic_category;
      typedef mpq_class field_extension_type;
    };

    template<> inline int convert_to<int>(const Numeric::Integer& n) { return n.get_si(); }
    template<> inline long convert_to<long>(const Numeric::Integer& n) { return n.get_si(); }
    
    template<> inline int conv_down<int>(const int& n) { return n; }
    template<> inline int conv_up<int>(const int& n) { return n; }
    

    //! \name %Integer arithmetical functions
    //@{
    //! \ingroup Numeric
    /*! \brief Factorial. */
    template<typename N>
    inline
    N 
    factorial(const N& n) {
      uint result=1;
      for(N i=1; i!=n; ++i) {
        result*=i;
      }
      return result;
    }


    /*! \brief Greatest common divisor. */
    template <typename N>
    inline 
    N 
    gcd(const N &a, const N &b)
    {
      N c=a%b;
      if (c==0) { 
        return b;
      }
      return (gcd(b, c));
    }
  
  
    /*! \brief Least common multiple. */
    template <typename N>
    inline 
    N 
    lcm(const N &a, const N &b) {
      return ((a*b)/gcd(a,b));
    }
    //@}

    //! \name %Integer bit-shift functions
    //@{
    /*! \brief The integer power \f$2^n\f$. */
    template<typename N>
    inline
    N 
    pow_two(const N& n) {
      return 1<<n;
    }

    /*! \brief The floor of the logarithm of \a n in base 2. */
    template<typename N>
    inline
    N 
    log_two_floor(const N& n) {
      assert(n>=1);
      N r=0;
      N y=n;
      while(y>=n) {
        y/=2;
        r+=1;
      }
      return r;
    }


    /*! \brief The ceiling of the logarithm of \a n in base 2. */
    template<typename N>
    inline
    N 
    log_two_ceil(const N& n) {
      assert(n>=1);
      N r=0;
      N y=n;
      while(y>1) {
        y/=2;
        r+=1;
      }
      return r;
    }

    //@}
  }
}

#endif /* _ARIADNE_INTEGER_H */
