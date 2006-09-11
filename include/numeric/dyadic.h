/***************************************************************************
 *            dyadic.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
 ****************************************************************************/

/*! \file dyadic.h
 *  \brief Type definitions and conversion operators for dyadic numbers.
 */

#ifndef _ARIADNE_DYADIC_H
#define _ARIADNE_DYADIC_H

#include "../declarations.h"
#include "../base/dyadic.h"
#include "../numeric/numerical_traits.h"

namespace Ariadne {
  namespace Numeric {

#ifdef DOXYGEN
    /*!\ingroup Numeric
     * \brief A dyadic rational (i.e. of form \f$m/2^n\f$).
     * 
     * A element of the ring of dyadic rationals.
     * Must allow denotation of any dyadic rational.
     * May be created without loss of precision from any integral or floating point type,
     * or from any rational of the form m/2^n.
     * May be created without loss of precision from any integral or floating point type,
     * or from any rational of the form m/2^n.
     *
     * Currently implemented using a modification of the Synaps dyadic class.
     *
     * FIXME: mpf_class does not implement addition, subtraction and multiplication exactly.
     */
    class Dyadic { };
#else
    typedef Synaps::dyadic Dyadic;
#endif

    template<> class numerical_traits<Dyadic> {
     public:
      typedef ring_tag algebraic_category;
      typedef Rational field_extension_type;
    };

    template<> inline std::string name<Numeric::Dyadic>() { return "Dyadic"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Dyadic> >() { return "Interval<Dyadic>"; }
 
    inline Dyadic mantissa(const Dyadic& num) {
      return num.mantissa();
    }
  
    inline int exponent(const Dyadic& num) {
      return num.exponent();
    }
  
    inline int precision(const Dyadic& num) {
      return num.precision();
    }
  
    inline Integer numerator(const Dyadic& num) {
      return num.numerator();
    }
  
    inline Integer denominator(const Dyadic& num) {
      return num.denominator();
    }
 
    
    

  } 
}

#endif // SYNAPS_DYADIC_H
