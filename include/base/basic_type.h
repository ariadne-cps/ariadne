/***************************************************************************
 *            basic_type.h
 *
 *  31 January 2006
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

/*! \file basic_type.h
 *  \brief Basic type definitions.
 */

#ifndef _ARIADNE_BASIC_TYPE_H
#define _ARIADNE_BASIC_TYPE_H

namespace Ariadne {
  namespace Base {
    template<typename R, size_t N=0> class array;
    template<typename R> class array_vector;
    template<typename R> class sequence;

    template<typename R> class Interval;
    class BinaryWord;
    class BinaryTree;
  }
}

namespace Ariadne {
  /*! \brief An unsigned integral type used to represent a coordinate in state space. */
  typedef unsigned short dimension_type;
  /*! \brief An unsigned integral type used to represent the size of a list. */
  typedef size_t size_type;
  /*! \brief An signed integral type used to represent the position in a list with positive and negative indices. */
  typedef int index_type;
  
  using Base::array;
  using Base::array_vector;
  using Base::sequence;

  using Base::Interval;
  using Base::BinaryWord;
  using Base::BinaryTree;

  /*! \brief An array of boolean values. */
  typedef Base::array<bool> BooleanArray;
  /*! \brief An array of unsigned integer values. */
  typedef Base::array<size_type> SizeArray;
  /*! \brief An array of integer values. */
  typedef Base::array<index_type> IndexArray;

}

#endif /* _ARIADNE_BASIC_TYPE_H */
