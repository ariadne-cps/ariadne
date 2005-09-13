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
 *  \brief Header for non-templated operations on grids.
 */

#ifndef GRID_OPERATIONS_H
#define GRID_OPERATIONS_H

#include <vector>
#include "array.h"
#include "rectangle.h"

namespace Ariadne {
   namespace Geometry {

    typedef unsigned short dimension_type;
    typedef size_t size_type;
    typedef int index_type;

    /*!\brief An array of bool, to be used as a mask. */
    typedef std::vector<bool> BooleanArray;

    /*!\brief An of integers representing an index in a grid. */
    typedef array<index_type> IndexArray;
     /*!\brief An of integers representing the size of a rectangle in a grid. */
   typedef array<size_type> SizeArray;
    /*!\brief An array of integers, representing a cell in a grid. */
    typedef array<index_type> IntegerCell;
    /*!\brief A rectangle in an integer grid. */
    typedef Rectangle<index_type> IntegerRectangle;
     /*!\brief A list of arrays of integers of the same size, representing cells in a grid. */
    typedef array_vector<index_type> IntegerCellList;
    /*!\brief A list of arrays of integers of the same size, representing rectangles in a grid. */
    typedef array_vector<index_type> IntegerRectangleList;

    inline size_type inner_product(const array<size_type>& a1, const array<size_type>& a2);

    inline BooleanArray& operator&=(BooleanArray& v1, const BooleanArray& v2);
    inline BooleanArray& operator|=(BooleanArray& v1, const BooleanArray& v2);
    inline BooleanArray& operator-=(BooleanArray& v1, const BooleanArray& v2);

    inline BooleanArray operator&(const BooleanArray& v1, const BooleanArray& v2);
    inline BooleanArray operator|(const BooleanArray& v1, const BooleanArray& v2);
    inline BooleanArray operator-(const BooleanArray& v1, const BooleanArray& v2);

    inline bool operator<(const IndexArray&, const IndexArray&);

    inline void compute_cell_list_bounds(IndexArray& l, IndexArray& u, IntegerCellList cl);
    inline void compute_rectangle_list_bounds(IndexArray& l, IndexArray& u, IntegerRectangleList cl);

    /*!\internal TODO:
     * unique sort cell list
     * subset remove rectangle list
     * compute redundant indices in rectangle list
     * compute redundant indices of mask array
     *
     * convert rectangle list to cell list
     * convert rectangle list to mask array
     * convert cell list to mask array
     * convert mask array to cell list
     *
     * transform indices in cell list
     * transform indices in rectangle list
     * transform indices in mask array
     *
     * compute bounding box of cell list DONE
     * compute bounding box of rectangle list DONE
     * compute bounding box of mask array
     *
     * sorted union cell list
     * union rectangle list
     * union mask array DONE
     *
     * intersection cell list
     * intersection rectangle list
     * intersection mask array DONE
     *
     * subtraction cell list
     * subtraction mask array DONE
     *
     * inplace subtraction cell list
     * inplace union mask array with all
     * inplace intersection mask array with all
     * inplace subtraction mask array with all
     *
     * mask array subset
     * mask array interior subset
     * mask array subset of rectangle
     */


  }
}

#endif /* _GRID_OPERATIONS_H */
