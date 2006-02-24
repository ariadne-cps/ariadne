/***************************************************************************
 *            unit_grid_set.h
 *
 *  Copyright  2005,6  Alberto Casagrande, Pieter Collins
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

/*! \file unit_grid_set.h
 *  \brief Sets on an integer grid.
 */

#ifndef _ARIADNE_UNIT_GRID_SET_H
#define _ARIADNE_UNIT_GRID_SET_H

#include <vector>
#include <iosfwd>

#include "base/array.h"
#include "base/basic_type.h"

#include "grid_operations.h"

namespace Ariadne {
  namespace Geometry {
    /*!\brief An array of bool, to be used as a mask. */
    typedef array<bool> BooleanArray;

    /*!\brief An of integers representing a cell in a grid. */
    typedef array<index_type> IndexArray;
     /*!\brief An of integers representing the size of a rectangle in a grid. */
    typedef array<size_type> SizeArray;

    /*!\brief An array of integers, representing a cell in a grid. */
    class UnitGridCell;
    /*!\brief An of integers representing an index in a grid. */
    class UnitGridRectangle;
    class UnitGridRectangleIterator;

    /*!\brief A list of arrays of integers of the same size, representing cells in a grid. */
    typedef array_vector<index_type> UnitCellList;
    /*!\brief A list of arrays of integers of the same size, representing rectangles in a grid. */
    typedef array_vector<index_type> UnitRectangleList;

    /*!\brief A list of arrays of integers of the same size, representing rectangles in a grid. */
    class UnitGridMaskSet;
    class UnitGridMaskSetIterator;

    UnitGridRectangle regular_intersection(const UnitGridRectangle&, const UnitGridRectangle&);
    bool subset(const UnitGridRectangle&, const UnitGridRectangle&);
    bool interiors_intersect(const UnitGridRectangle&, const UnitGridRectangle&);

    bool interiors_intersect(const UnitGridRectangle&, const UnitGridMaskSet&);
    bool subset(const UnitGridRectangle&, const UnitGridMaskSet&);
    
    std::ostream& operator<<(std::ostream& os, const UnitGridRectangle&);
    

    /*! \brief A cell in a unit grid. */
    class UnitGridCell {
      friend class UnitGridRectangleIterator;
      friend class UnitGridMaskSet;
     public:
      UnitGridCell(dimension_type n) : _lower(n) { }
      UnitGridCell(const IndexArray& l) : _lower(l) { }
      UnitGridCell(const UnitGridCell& c) : _lower(c._lower) { }
      
      bool operator==(const UnitGridCell& c) const { 
        return this->_lower==c._lower; }
        
      dimension_type dimension() const { return this->_lower.size(); }
      Interval<index_type> operator[](dimension_type i) {
        return Interval<index_type>(this->_lower[i],this->_lower[i]+1); } 
      index_type lower_bound(dimension_type i) const { return this->_lower[i]; }
      index_type upper_bound(dimension_type i) const { return this->_lower[i]+1; }

      const IndexArray& position() const { return this->_lower; }
      IndexArray lower() const { return this->_lower; }
      IndexArray upper() const { 
        IndexArray result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->upper_bound(i); }
        return result;
      }
     private:
      IndexArray _lower;
    };
    
    /*! \brief A block of indicejos in a grid. */
    class UnitGridRectangle {
      friend class UnitGridRectangleIterator;
      friend class UnitGridMaskSet;
     public:
      typedef UnitGridRectangleIterator const_iterator;
      UnitGridRectangle(dimension_type n) : _lower(n), _upper(n) { }
      UnitGridRectangle(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u) { }
      UnitGridRectangle(const UnitGridCell& c)
        : _lower(c.lower()), _upper(c.upper()) { }
      UnitGridRectangle(const UnitGridRectangle& r)
        : _lower(r._lower), _upper(r._upper) { }
      
      bool operator==(const UnitGridRectangle& other) const { 
        return this->_lower==other._lower && this->_upper==other._upper; }
        
      dimension_type dimension() const { return this->_lower.size(); }
      bool empty() const;
      
      Interval<index_type> operator[](dimension_type i) const {
        return Interval<index_type>(this->_lower[i],this->_upper[i]); } 
      index_type lower_bound(dimension_type i) const { return this->_lower[i]; }
      index_type upper_bound(dimension_type i) const { return this->_upper[i]; }

      void set_lower_bound(dimension_type i, index_type n) { this->_lower[i]=n; }
      void set_upper_bound(dimension_type i, index_type n) { this->_upper[i]=n; }

      const IndexArray& lower() const { return this->_lower; }
      const IndexArray& upper() const { return this->_upper; }
      SizeArray sizes() const { return this->_upper-this->_lower; }
      SizeArray strides() const { return compute_strides(this->sizes()); }
      size_type size() const { return compute_strides(this->sizes())[dimension()]; }
      
      const_iterator begin() const;
      const_iterator end() const;
     private:
      IndexArray _lower;
      IndexArray _upper;
    };
    
    /*! \brief An iterator for positions in rectangular piece of a grid. */
    class UnitGridRectangleIterator {
     public:
      UnitGridRectangleIterator(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u), _position(l) { }
      UnitGridRectangleIterator(const IndexArray& l, const IndexArray& u, const IndexArray& p)
        : _lower(l), _upper(u), _position(p) { }
      UnitGridCell operator*() const { return this->dereference(); }
      UnitGridRectangleIterator& operator++() { this->increment(); return *this; }
      
      bool operator==(const UnitGridRectangleIterator& other) const {
        return this->equal(other); }
      bool operator!=(const UnitGridRectangleIterator& other) const {
        return !this->equal(other); }
     private:
      bool equal(const UnitGridRectangleIterator& other) const {
        return (this->_position==other._position) 
          && (this->_lower==other._lower) && (this->_upper==other._upper);
      }
      void increment() {
        dimension_type d=0;
        _position[d]+=1;
        while(_position[d]==_upper[d] && (d+1u)!=_position.size() ) {
          _position[d]=_lower[d];
          d+=1;
          _position[d]+=1;
        }
      }
      UnitGridCell dereference() const { return UnitGridCell(_position); }
     private:
      dimension_type dimension() const { return _position.size(); }
      const IndexArray& _lower;
      const IndexArray& _upper;
      IndexArray _position;
    };

    inline
    UnitGridRectangle::const_iterator 
    UnitGridRectangle::begin() const 
    { 
      return const_iterator(this->_lower, this->_upper, this->_lower);
    }
    
    inline
    UnitGridRectangle::const_iterator 
    UnitGridRectangle::end() const 
    { 
      IndexArray end_position=this->_lower;
      end_position[this->dimension()-1]=_upper[this->dimension()-1];
      return const_iterator(this->_lower, this->_upper, end_position);
    }
    
     /*!\brief A list of arrays of integers of the same size, representing rectangles in a grid. */
    class UnitGridMaskSet {
     public:
      typedef UnitGridMaskSetIterator iterator;
      typedef UnitGridMaskSetIterator const_iterator;
     public:
      UnitGridMaskSet(const UnitGridRectangle& bb) 
        : _bounds(bb), _mask(bb.size(),false) { this->_compute_cached_attributes(); }
      UnitGridMaskSet(const UnitGridRectangle& bb, const BooleanArray& ma) 
        : _bounds(bb), _mask(ma) { this->_compute_cached_attributes(); }

      UnitGridMaskSet(const UnitGridMaskSet& ms);
      
      const UnitGridRectangle& bounds() const { return _bounds; }
      const BooleanArray& mask() const { return _mask; };
      dimension_type dimension() const { return _bounds.dimension(); }
      size_type capacity() const { return _mask.size(); }
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }
      size_type empty() const { return this->size()==0; }
      const IndexArray& lower() const { return _lower; }
      const IndexArray& upper() const { return _upper; }
      const SizeArray& sizes() const { return _sizes; }
      const SizeArray& strides() const { return _strides; }
      
      void adjoin(const UnitGridCell& c) { 
        compute_cell_mask(&_mask,_strides,_lower,c.position()); 
      }
      void adjoin(const UnitGridRectangle& r) { 
        compute_rectangle_mask(&_mask,_strides,_lower,
                               r.lower(),r.upper());
      }
      /*! \brief Adjoins a UnitGridMaskSet to the set. */
      void adjoin(const UnitGridMaskSet& ms) {
        assert(ms._bounds==this->_bounds);
        _mask |= ms._mask;
      }

      /*! \brief Adjoins a GridCellListSet to the set. */
      void adjoin_cells(const UnitCellList& ucl) {
        compute_cell_list_mask(&_mask,_strides,_lower,ucl);
      }

      /*! \brief Adjoins a GridRectangleListSet to the set. */
      void adjoin_rectangles(const UnitRectangleList& url) {
        compute_rectangle_list_mask(&_mask,_strides,_lower,url);
      }
      
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const;
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const;
     private:
      void _compute_cached_attributes();
     private:
      UnitGridRectangle _bounds;
      IndexArray _lower;
      IndexArray _upper;
      SizeArray _sizes;
      SizeArray _strides;
      BooleanArray _mask;
    };

        
     
    //TODO: Replace by a general-purpose mask iterator
    class UnitGridMaskSetIterator {
      typedef UnitGridCell value_type;
      typedef UnitGridCell reference;
     public:
      UnitGridMaskSetIterator(UnitGridRectangle::const_iterator ri, 
                              BooleanArray::const_iterator mi, 
                              BooleanArray::const_iterator me) 
        : _cell_iter(ri), _mask_iter(mi), _mask_end(me) 
      { 
        while(_mask_iter!=_mask_end && !*_mask_iter) { 
          ++_mask_iter; ++_cell_iter; 
        }
      }
      
      bool operator==(const UnitGridMaskSetIterator& other) const { 
        return this->equal(other); }
      bool operator!=(const UnitGridMaskSetIterator& other) const { 
        return !this->equal(other); }
      UnitGridCell operator*() const { 
        return this->dereference(); }
      UnitGridMaskSetIterator& operator++() { 
        this->increment(); return *this; }
     private:
      bool equal(const UnitGridMaskSetIterator& other) const {
        return this->_cell_iter==other._cell_iter && this->_mask_iter==other._mask_iter;
      }
      UnitGridCell dereference() const { 
        return _cell_iter.operator*(); 
      }
      void increment() { 
        do { 
          ++_mask_iter; ++_cell_iter; } 
        while(_mask_iter!=_mask_end && !*_mask_iter);
      }
     private:
      UnitGridRectangle::const_iterator _cell_iter;
      BooleanArray::const_iterator _mask_iter;
      BooleanArray::const_iterator _mask_end;
    };

    inline
    UnitGridMaskSet::const_iterator 
    UnitGridMaskSet::begin() const 
    { 
      return const_iterator(_bounds.begin(),_mask.begin(),_mask.end()); 
    }
    
    inline
    UnitGridMaskSet::const_iterator 
    UnitGridMaskSet::end() const 
    { 
      return const_iterator(_bounds.end(),_mask.end(),_mask.end());
    }

  }
}

#endif /* _ARIADNE_UNIT_GRID_SET_H */
