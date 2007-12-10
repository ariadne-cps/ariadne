/***************************************************************************
 *            numeric/floatmp.class.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file numeric/floatmp.class.h
 *  \brief Definition of multiple-precision floating point class.
 */

#ifndef ARIADNE_NUMERIC_FLOATMP_CLASS_H
#define ARIADNE_NUMERIC_FLOATMP_CLASS_H

#include <mpfr.h>
#include <gmp.h>

#include "numeric/expression.h"

namespace Ariadne {
  namespace Numeric {

    class mpfr;
    
    template<class T> class Float;
    template<class R> class Interval;

    typedef Float<mpfr> FloatMP;
    typedef Interval<FloatMP> IntervalMP;



    /*!\ingroup Numeric
     * \brief A multiple-precision floating-point type.
     *  
     * Currently implemented using mpfr_t from the MPFR library.
     */
    template<> 
    class Float<mpfr>
      : public Value< Float<mpfr> >
    {
     public:
      ~Float();
      Float();
      Float(const int& n);
      Float(const uint& n);
      Float(const double& x);
      Float(const Integer& n);
      Float(const FloatMP& x);

      explicit Float(const std::string& x);

      FloatMP& operator=(const int& n);
      FloatMP& operator=(const uint& n);
      FloatMP& operator=(const double& x);
      FloatMP& operator=(const Integer& n);
      FloatMP& operator=(const FloatMP& x);

      template<class E> Float(const Expression<E>& e);
      template<class E> FloatMP& operator=(const Expression<E>& e);

      template<class X, class Rnd> Float(const X& x, Rnd rnd);
      template<class E, class Rnd> Float(const Expression<E>& e, Rnd rnd);

      static uint default_precision();
      static void set_default_precision(uint p);
      void set_precision(uint p);
      uint precision() const;
     public:
      mpfr_t _value;
    };
    
    std::ostream& operator<<(std::ostream& os, const FloatMP& x);
    std::istream& operator>>(std::istream& is, FloatMP& x);
      
  }
}

#endif /* ARIADNE_NUMERIC_FLOATMP_CLASS_H */
