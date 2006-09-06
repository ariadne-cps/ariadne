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
 *  \brief Type definitions and conversion operators for fundamental Ariadne types.
 */

#ifndef _ARIADNE_INTEGER_H
#define _ARIADNE_INTEGER_H

#include <gmpxx.h>

#include "../declarations.h"
#include "../numeric/numerical_traits.h"

namespace Ariadne {
  namespace Numeric {
#ifdef DOXYGEN
    /*! \brief An integer.
     *  \ingroup Numeric
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
    
  }
}

#endif /* _ARIADNE_INTEGER_H */
