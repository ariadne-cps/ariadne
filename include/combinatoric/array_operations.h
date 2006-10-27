/***************************************************************************
 *            array_operations.h
 *
 *  22 June 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

/*! \file array_operations.h
 *  \brief Operations on arrays.
 */

#ifndef _ARIADNE_ARRAY_OPERATIONS_H
#define _ARIADNE_ARRAY_OPERATIONS_H

#include "../declarations.h"


namespace Ariadne {
  namespace Base {
    /*! \brief Inner product. */
    size_type inner_product(const SizeArray& a1, const SizeArray& a2);
    /*! \brief Inner product. */
    index_type inner_product(const IndexArray& a1, const IndexArray& a2);

    /*! \brief Compare two arrays using the lexicographic total ordering. */
    bool lexicographic_less(const IndexArray&, const IndexArray&);
    /*! \brief Compare two arrays using the lexicographic total ordering. */
    bool lexicographic_less_equal(const IndexArray&, const IndexArray&);
    /*! \brief Compare two arrays using the componentwise partial ordering. */
    bool coordinate_less_equal(const IndexArray&, const IndexArray&);
  
    /*! \brief Assigns the max of a and b to a. */
    void assign_max(IndexArray& a, const IndexArray& l);

    /*! \brief Assigns the minimum of a and b componentwise to a. */
    void assign_min(IndexArray& a, const IndexArray& u);


    /*! \brief Compute the sum of an index array and a size. */
    IndexArray operator+(const IndexArray& l, const SizeArray& s);

    /*! \brief Compute a positive offset from two index sets. */
    IndexArray operator-(const IndexArray& u, const IndexArray& l);

    /*! \brief Inplace boolean and. */
    BooleanArray& operator&=(BooleanArray& v1, const BooleanArray& v2);
    /*! \brief Inplace boolean or. */
    BooleanArray& operator|=(BooleanArray& v1, const BooleanArray& v2);
    /*! \brief Inplace boolean minus. */
    BooleanArray& operator-=(BooleanArray& v1, const BooleanArray& v2);
    
    /*! \brief Boolean and. */
    BooleanArray operator&(const BooleanArray& v1, const BooleanArray& v2);
    /*! \brief Boolean or. */
    BooleanArray operator|(const BooleanArray& v1, const BooleanArray& v2);
    /*! \brief Boolean minus. */
    BooleanArray operator-(const BooleanArray& v1, const BooleanArray& v2);

    /*! \brief True if v1-v2 is all zeros. */
    bool operator<=(const BooleanArray& v1, const BooleanArray& v2);
    

  }

}


#endif /* _ARIADNE_ARRAY_OPERATIONS_H */
