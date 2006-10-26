/***************************************************************************
 *            real_typedef.h
 *
 *  06 Feb 2006
 *  Copyright  2005  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file real_typedef.h
 *  \brief Defines the real number type used by the library and the Python interface.
 */

#ifndef _ARIADNE_REAL_TYPEDEF_H
#define _ARIADNE_REAL_TYPEDEF_H


#ifndef DOUBLE_REAL 
#define MPFLOAT_REAL
#endif

#ifdef DOUBLE_REAL
#define REAL_IS_A_FIELD
#include "numeric/float64.h"
namespace Ariadne {
  typedef Numeric::Float64 Real;
}
#endif

#ifdef MPFLOAT_REAL
#undef REAL_IS_A_FIELD
#include "numeric/mpfloat.h"
namespace Ariadne {
  typedef Numeric::MPFloat Real;
}
#endif

namespace Ariadne {
  typedef Ariadne::Numeric::traits<Real>::arithmetic_type Field;
}

#endif /* _ARIADNE_REAL_TYPEDEF_H */
