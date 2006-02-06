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
#include "array.h"
#include "basic_type.h"
#include "rectangle.h"

namespace Ariadne {
  namespace Geometry {
    /*!\brief An array of bool, to be used as a mask. */
    typedef array<bool> BooleanArray;

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

    size_type inner_product(const array<size_type>& a1, const array<size_type>& a2);

    BooleanArray& operator&=(BooleanArray& v1, const BooleanArray& v2);
    BooleanArray& operator|=(BooleanArray& v1, const BooleanArray& v2);
    BooleanArray& operator-=(BooleanArray& v1, const BooleanArray& v2);
    
    BooleanArray operator&(const BooleanArray& v1, const BooleanArray& v2);
    BooleanArray operator|(const BooleanArray& v1, const BooleanArray& v2);
    BooleanArray operator-(const BooleanArray& v1, const BooleanArray& v2);
    
    bool operator<(const IndexArray&, const IndexArray&);
    
    /*! Compute the sum of an index array and a size. */
    IndexArray operator+(const IndexArray& l, const SizeArray& s);

    /*! Compute a positive offset from two index sets */
    SizeArray operator-(const IndexArray& u, const IndexArray& l);
   
    void compute_rectangle_list_bounds(IndexArray& l, IndexArray& u, IntegerRectangleList cl);

    /*! Compute the index of a position in a grid. */
    size_type compute_index(const IndexArray& pos, const IndexArray& lower, const SizeArray& strides);

    /*! Compute the position of an index in a grid. */
    IndexArray compute_position(size_type index, const IndexArray& lower, const SizeArray& strides);

    /*! Compute upper and lower bounds of the cell list cl. */
      void compute_cell_list_bounds(IndexArray* lptr, IndexArray* uptr, IntegerCellList cl);

    /*! Compute strides from a list of sizes. */
    SizeArray compute_strides(const SizeArray& s);

    /*! Compute lower bounds of the cell list cl. */
    void compute_rectangle_list_lower_bound(IndexArray* lptr, const IntegerRectangleList& rl);

    /* Compute upper bounds of the cell list cl. */
    void compute_rectangle_list_upper_bound(IndexArray* uptr, const IntegerRectangleList& rl);

    /*! Compute upper and lower bounds of the cell list cl. */
    void compute_rectangle_list_bounds(IndexArray* lptr, IndexArray* uptr, const IntegerRectangleList& rl);

    void append_to_cell_list(IntegerCellList* clptr, const IndexArray& lower, const SizeArray& strides, const BooleanArray& mask);
    void append_to_cell_list(IntegerCellList* clptr, const IndexArray& lower, const IndexArray& upper);
    void append_to_cell_list(IntegerCellList* clptr, const IntegerRectangleList rl);

    void compute_cell_mask(BooleanArray* maptr, const SizeArray& grid_strides, const IndexArray& grid_lower, const IndexArray& position);
    void compute_cell_list_mask(BooleanArray* maptr, const SizeArray& grid_strides, const IndexArray& grid_lower, const IntegerCellList& cl);
    void compute_rectangle_mask(BooleanArray* maptr, const SizeArray& grid_strides, const IndexArray& grid_lower, 
                                const IndexArray& lower, const IndexArray& upper);
    void compute_rectangle_list_mask(BooleanArray* maptr, const SizeArray& grid_strides, 
                                     const IndexArray& grid_lower, const IntegerRectangleList& rl);

    void translate_rectangle_coordinates(IntegerRectangleList* torlptr, const IntegerRectangleList& frrl, array< std::vector<index_type> > tr);
    void translate_cell_coordinates(IntegerRectangleList* torlptr, const IntegerCellList& frcl, array< std::vector<index_type> > tr);

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

#endif /* _ARIADNE_GRID_OPERATIONS_H */
