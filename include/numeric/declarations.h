/***************************************************************************
 *            numeric/declarations.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file numeric/declarations.h
 *  \brief Forward declarations of classes in the Numeric module.
 */

#ifndef ARIADNE_NUMERIC_DECLARATIONS_H
#define ARIADNE_NUMERIC_DECLARATIONS_H

#include "rounding.h"

namespace Ariadne { 
  

    class mpfr;

    class Integer;
    class Rational;
    template<class T> class Float;
    template<class T> class ApproximateFloat;
    template<class T> class ErrorFloat;
    template<class R> class Interval;

    typedef Float<double> Float64;
    typedef Float<mpfr> FloatMP;
  
  
} // namespace Ariadne


#endif /* ARIADNE_NUMERIC_DECLARATIONS_H */
