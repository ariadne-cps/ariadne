/***************************************************************************
 *            builtin.h
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

/*! \file builtin.h
 *  \brief Built-in and standard classes and types.
 */

#ifndef ARIADNE_BUILTIN_H
#define ARIADNE_BUILTIN_H

#ifdef DOXYGEN

//! \name Built-in types
//@{
/*!\brief Built-in logical type. */
class bool { };
/*!\brief Built-in integer type. */
class int { };
/*!\brief Built-in unsigned integer type. */
class uint { };
/*!\brief Built-in double-precision floating-point type. */
class double { };
/*!\brief Built-in character type. */
class char { };
//@}

namespace std {
  //! \name Standard library classes
  //@{
  /*!\brief Standard %string type. */
  class string { };
  //@}
}

#endif



#endif /* ARIADNE_BUILTIN_H */
