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
    /*!\brief An array of bool, to be used as a mask. */
    typedef array<bool> BooleanArray;

    /*!\brief An of integers representing an index in a grid. */
    typedef array<index_type> IndexArray;
     /*!\brief An of integers representing the size of a rectangle in a grid. */
    typedef array<size_type> SizeArray;
    /*!\brief An array of integers, representing a cell in a grid. */
    typedef array<index_type> IntegerCell;
    /*!\brief An of integers representing an index in a grid. */
    class IndexBlock;
     /*!\brief A list of arrays of integers of the same size, representing cells in a grid. */
    typedef array_vector<index_type> IndexArrayList;
    /*!\brief A list of arrays of integers of the same size, representing rectangles in a grid. */
    typedef array_vector<index_type> IndexBlockList;

    size_type inner_product(const array<size_type>& a1, const array<size_type>& a2);

    BooleanArray& operator&=(BooleanArray& v1, const BooleanArray& v2);
    BooleanArray& operator|=(BooleanArray& v1, const BooleanArray& v2);
    BooleanArray& operator-=(BooleanArray& v1, const BooleanArray& v2);
    
    BooleanArray operator&(const BooleanArray& v1, const BooleanArray& v2);
    BooleanArray operator|(const BooleanArray& v1, const BooleanArray& v2);
    BooleanArray operator-(const BooleanArray& v1, const BooleanArray& v2);
    
    bool operator<(const IndexArray&, const IndexArray&);
    
    /*! Returns true if v1-v2 is all zeros. */
    bool operator<=(const BooleanArray& v1, const BooleanArray& v2);
    
    /*! Compute the sum of an index array and a size. */
    IndexArray operator+(const IndexArray& l, const SizeArray& s);

    /*! Compute a positive offset from two index sets */
    SizeArray operator-(const IndexArray& u, const IndexArray& l);
   
    /*! Compute the index of a position in a grid. */
    size_type compute_index(const IndexArray& pos, const IndexArray& lower, const SizeArray& strides);

    /*! Compute the position of an index in a grid. */
    IndexArray compute_position(size_type index, const IndexArray& lower, const SizeArray& strides);

    /*! Compute upper and lower bounds of the cell list cl. */
    IndexBlock compute_cell_list_bounds(const IndexArrayList& cl);

    /*! Compute strides from a list of sizes. */
    SizeArray compute_strides(const SizeArray& s);

    /*! Compute upper and lower bounds of the rectangle list rl. */
    IndexBlock compute_rectangle_list_bounds(const IndexBlockList& rl);

    void append_to_cell_list(IndexArrayList* clptr, const IndexArray& lower, const SizeArray& strides, const BooleanArray& mask);
    void append_to_cell_list(IndexArrayList* clptr, const IndexArray& lower, const IndexArray& upper);
    void append_to_cell_list(IndexArrayList* clptr, const IndexBlockList rl);

    void compute_cell_mask(BooleanArray* maptr, const SizeArray& grid_strides, const IndexArray& grid_lower, const IndexArray& position);
    void compute_cell_list_mask(BooleanArray* maptr, const SizeArray& grid_strides, const IndexArray& grid_lower, const IndexArrayList& cl);
    void compute_rectangle_mask(BooleanArray* maptr, const SizeArray& grid_strides, const IndexArray& grid_lower, 
                                const IndexArray& lower, const IndexArray& upper);
    void compute_rectangle_list_mask(BooleanArray* maptr, const SizeArray& grid_strides, 
                                     const IndexArray& grid_lower, const IndexBlockList& rl);

    void translate_rectangle_coordinates(IndexBlockList* torlptr, const IndexBlockList& frrl, array< std::vector<index_type> > tr);
    void translate_cell_coordinates(IndexBlockList* torlptr, const IndexArrayList& frcl, array< std::vector<index_type> > tr);



    /*! \brief An iterator for positions in rectangular piece of a grid. */
    class GridPositionIterator {
     public:
      GridPositionIterator(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u), _position(l) { }
      GridPositionIterator(const IndexArray& l, const IndexArray& u, const IndexArray& p)
        : _lower(l), _upper(u), _position(p) { }
      const IndexArray& operator*() const { return _position; }
      GridPositionIterator& operator++() {
        dimension_type d=0;
        _position[d]+=1;
        while(_position[d]==_upper[d] && (d+1u)!=_position.size() ) {
          _position[d]=_lower[d];
          d+=1;
          _position[d]+=1;
        }
        return *this;
      }
      bool operator==(const GridPositionIterator& other) const {
        return (this->_position==other._position) 
          && (this->_lower==other._lower) && (this->_upper==other._upper);
      }
      bool operator!=(const GridPositionIterator& other) const {
        return !(*this==other);
      }
      bool end() const { return _position[dimension()-1]==_upper[dimension()-1]; }
     private:
      dimension_type dimension() const { return _position.size(); }
      const IndexArray& _lower;
      const IndexArray& _upper;
      IndexArray _position;
    };


    /*! \brief A block of indicejos in a grid. */
    class IndexBlock {
     public:
      typedef GridPositionIterator const_iterator;
      IndexBlock(dimension_type n) : _lower(n), _upper(n) { }
      IndexBlock(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u) { }
      IndexBlock(const IndexBlock& b)
        : _lower(b._lower), _upper(b._upper) { }
      
      bool operator==(const IndexBlock& b) { 
        return this->_lower==b._lower && this->_upper==b._upper; }
        
      dimension_type dimension() const { return _lower.size(); }
      index_type lower_bound(dimension_type i) const { return _lower[i]; }
      index_type upper_bound(dimension_type i) const { return _upper[i]; }

      void set_lower_bound(dimension_type i, index_type n) { _lower[i]=n; }
      void set_upper_bound(dimension_type i, index_type n) { _upper[i]=n; }

      const IndexArray& lower() const { return _lower; }
      const IndexArray& upper() const { return _upper; }
      SizeArray sizes() const { return _upper-_lower; }
      SizeArray strides() const { return compute_strides(sizes()); }
      
      GridPositionIterator begin() const { 
        return GridPositionIterator(_lower, _upper,_lower);
      }
      GridPositionIterator end() const { 
        IndexArray end_position=_lower;
        end_position[this->dimension()-1]=_upper[this->dimension()-1];
        return GridPositionIterator(_lower, _upper,end_position);
      }
     private:
      IndexArray _lower;
      IndexArray _upper;
    };
    
    std::ostream&
    operator<<(std::ostream&, const IndexBlock&);
    
    
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
