/***************************************************************************
 *            types.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file types.h
 *  \brief Basic Ariadne typedef's
 */

#ifndef ARIADNE_TYPES_H
#define ARIADNE_TYPES_H

#include <cstdlib>

#include "array.decl.h"

namespace Ariadne { 

  namespace Numeric {
    class Integer;
    class Rational;
  }
  
  
  namespace Base {
    //! \name Basic types
    //@{
    //! \ingroup Base
    /*! \brief The type used for a unique identifier or key. */
    typedef size_t id_type;
    /*! \brief The type of a machine byte. */
    typedef unsigned char byte_type; 
    /*! \brief An unsigned integral type used to represent the size of a list. */
    typedef size_t size_type;
    /*! \brief An signed integral type used to represent the position in a list with positive and negative indices. */
    typedef int index_type;
    /*! \brief An unsigned integral type used to represent a coordinate in state space. */
    typedef unsigned short dimension_type;
    /*! \brief The type used to describe evolution time for discrete time systems. */
    typedef Numeric::Integer discrete_time_type;
    /*! \brief The type used to describe evolution time. */
    typedef Numeric::Rational time_type;
    //@}
  
    //! \name Array types
    //@{
    //! \ingroup Array
    /*! \brief An array of boolean values. */
    typedef array<bool> BooleanArray;
    /*! \brief An array of unsigned integer values. */
    typedef array<size_type> SizeArray;
    /*! \brief An array of integer values. */
    typedef array<index_type> IndexArray;
    //@}
  }
  
  using namespace Base;
  
}

#endif /* ARIADNE_TYPES_H */
