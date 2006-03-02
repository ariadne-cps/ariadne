/***************************************************************************
 *            grid_operations.h
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

/*! \file grid_operations.h
 *  \brief Non-templated operations on grids.
 */

#ifndef _ARIADNE_GRID_OPERATIONS_H
#define _ARIADNE_GRID_OPERATIONS_H

#include <vector>
#include <iosfwd>

#include "base/array.h"
#include "base/basic_type.h"

#include "geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {
    size_type inner_product(const array<size_type>& a1, const array<size_type>& a2);

    BooleanArray& operator&=(BooleanArray& v1, const BooleanArray& v2);
    BooleanArray& operator|=(BooleanArray& v1, const BooleanArray& v2);
    BooleanArray& operator-=(BooleanArray& v1, const BooleanArray& v2);
    
    BooleanArray operator&(const BooleanArray& v1, const BooleanArray& v2);
    BooleanArray operator|(const BooleanArray& v1, const BooleanArray& v2);
    BooleanArray operator-(const BooleanArray& v1, const BooleanArray& v2);
    
    bool lexicographic_order(const IndexArray&, const IndexArray&);
    bool coordinate_order(const IndexArray&, const IndexArray&);

    /*! Returns true if v1-v2 is all zeros. */
    bool operator<=(const BooleanArray& v1, const BooleanArray& v2);
    
    /*! Compute the sum of an index array and a size. */
    IndexArray operator+(const IndexArray& l, const SizeArray& s);

    /*! Compute a positive offset from two index sets */
    SizeArray operator-(const IndexArray& u, const IndexArray& l);
   
    /*! Assigns the max of a and b to a. */
    void assign_max(IndexArray& a, const IndexArray& l);

    /*! Assigns the minimum of a and b componentwise to a. */
    void assign_min(IndexArray& a, const IndexArray& u);

    /*! Compute strides from a list of sizes. */
    SizeArray compute_strides(const SizeArray& s);

    /*! Compute the index of a position in a grid. */
    size_type compute_index(const IndexArray& pos, const IndexArray& lower, const SizeArray& strides);

    /*! Compute the index of a position in a grid. */
    size_type compute_index(const IndexArray& pos, const SizeArray& strides, const index_type offset);

    /*! Compute the position of an index in a grid. */
    IndexArray compute_position(size_type index, const IndexArray& lower, const SizeArray& strides);


  }
}

#endif /* _ARIADNE_GRID_OPERATIONS_H */
