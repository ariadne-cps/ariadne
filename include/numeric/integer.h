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
#include <iostream>

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
  
    template<> inline int convert_to<int>(const Numeric::Integer& n) { return n.get_si(); }
    template<> inline long convert_to<long>(const Numeric::Integer& n) { return n.get_si(); }
    
    template<> inline int conv_down<int>(const int& n) { return n; }
    template<> inline int conv_up<int>(const int& n) { return n; }
    

    //! \name %Integer arithmetic
    //@{
    //! \ingroup Numeric
    template<> inline Integer min(const Integer& n1, const Integer& n2) {
      //std::cerr << "min<" << name<R>() << ">" << std::endl;
      return n1<=n2 ? n1 : n2;
    }
  
    template<> inline Integer max(const Integer& n1, const Integer& n2) {
      //std::cerr << "min<" << name<R>() << ">" << std::endl;
      return n1>=n2 ? n1 : n2;
    }
     
    /*! \brief Absolute value. */
    template<> inline 
    Integer abs(const Integer& n) {
      return n>=0 ? n : static_cast<Integer>(-n);
    }
  
    /*! \brief Unary negation. */
    template<> inline
    Integer neg(const Integer& n) {
      return -n;
    }
    
    /*! \brief Addition. */
    template<> inline
    Integer add(const Integer& n1,const Integer& n2) {
      return n1+n2;
    }

    
    /*! \brief Subtraction. */
    template<> inline
    Integer sub(const Integer& n1,const Integer& n2) {
      return n1-n2;
    }

    
    /*! \brief Multiplication. */
    template<> inline
    Integer mul(const Integer& n1,const Integer& n2) {
      return n1*n2;
    }
    
    /*! \brief The power of a real number type by an integer. */
    template<> inline 
    Integer pow(const Integer& n, const uint& i) {
      Integer r=1; Integer p=n; uint e=1;
      while(e<i) { if(e&i) { r*=p; } p*=p; e*=2; }
      return r; 
    }
    //@}
    
    //! \name %Integer functions
    //@{
    /*! \brief Factorial. */
    template<class N1, class N2> inline
    N1 factorial(const N2& n) {
      N1 result=1;
      for(N2 i=1; i!=n; ++i) { result*=i; }
      return result*n;
    }

    /*! \brief Factorial. */
    template<class N> inline
    N factorial(const N& n) {
      N result=1;
      if(n<=0) { return 1; }
      for(N i=1; i!=n; ++i) { result*=i; }
      return result*n;
    }

    /*! \brief Greatest common divisor. */
    template<class N> inline 
    N gcd(const N &a, const N &b) {
      N aa=a; N bb=b; N cc=aa%bb;
      while(cc!=0) { aa=bb; bb=cc; cc=aa%bb; }
      return bb;
    }

    /*! \brief Least common multiple. */
    template<class N> inline 
    N lcm(const N &a, const N &b) {
      return ((a*b)/gcd(a,b));
    }
    //@}

    //! \name %Integer bit-shift functions
    //@{
    /*! \brief The integer power \f$2^n\f$. */
    template<class N> inline
    N exp2(const N& n) {
      return 1<<n;
    }

    /*! \brief The floor of the logarithm of \a n in base 2. */
    template<class N> inline
    N log2_floor(const N& n) {
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
    template<class N> inline
    N log2_ceil(const N& n) {
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
