/***************************************************************************
 *            grid.h
 *
 *  18 January 2005
 *  Copyright  2004,2005  Alberto Casagrande, Pieter Collins
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
 
/*! \file grid.h
 *  \brief Cuboidal grids.
 */

#ifndef _GRID_H
#define _GRID_H

#include "interval.h"
#include "binary_word.h"
#include "array.h"
#include "sequence.h"

namespace Ariadne {	
  namespace Geometry {	
    template<typename R> class Rectangle;

    typedef unsigned short dimension_type;
    typedef size_t size_type;

    /*!\brief A finite, nonuniform grid of rectangles in Euclidean space.
     */
    template<typename R> 
    class FiniteGrid {
     public:
      FiniteGrid(const array< std::vector<R> >& sp)
	: _subdivision_positions(sp), _strides(sp.size()+1) { create(); }

      dimension_type dimension() const { return _subdivision_positions.size(); }
      size_type size() const { return _strides[dimension()]; }
      size_type size(dimension_type i) const { return _subdivision_positions.size()-1; }

      const array< std::vector<R> >& subdivision_positions() const { return _subdivision_positions; }
      const std::vector<R>& subdivision_positions(dimension_type d) const { return _subdivision_positions[d]; }
      R subdivision_position(dimension_type d, size_type n) const { return _subdivision_positions[d][n]; }
      
      inline size_type index(const array<size_type>& pos) const;
     private:
      void create();
     private:
      array< std::vector<R> > _subdivision_positions;
      array< size_type > _strides;
    };



    /*!\brief An infinite, uniform grid of rectangles in Euclidean space.
     */
    template<typename R> class InfiniteGrid {
     public:
      dimension_type dimension() const { return _subdivision_lengths.size(); }
      const array<R>& subdivision_lengths() const { return _subdivision_lengths; }
      R subdivision_length(dimension_type d) const { return _subdivision_lengths[d]; }
      R subdivision_position(dimension_type d, size_type n) const { return _subdivision_lengths[d] * n; }
     private:
      array<R> _subdivision_lengths;
    };



    /*!\brief A grid of rectangles in Euclidean space.
     */
    template<typename R> class PartitionGrid {
     public:
      PartitionGrid(const Rectangle<R>& bb, const sequence<dimension_type>& sc) 
	: _bounding_box(bb), _subdivision_coordinates(sc) { }

      dimension_type dimension() const { return _bounding_box.dimension(); }

      const Rectangle<R>& bounding_box() const { return _bounding_box; }
      const sequence<dimension_type>& subdivision_coordinates() const { return _subdivision_coordinates; }	
      dimension_type subdivision_coordinate(size_type n) const { return _subdivision_coordinates[n]; }	
     private:
      Rectangle<R> _bounding_box;
      sequence<dimension_type> _subdivision_coordinates;
    };


    template<typename R>
    void
    FiniteGrid<R>::create() 
    {
      _strides[0]=1;
      for(dimension_type i=0; i!=dimension(); ++i) {
	std::vector<R>& pos=_subdivision_positions[i];
	std::sort(pos.begin(),pos.end());
	typename std::vector<R>::iterator newend=std::unique(pos.begin(),pos.end());
	pos.resize(std::distance(pos.begin(),newend));
	_strides[i+1] = _strides[i] * size(i);
      }
    }

    template<typename R>
    inline 
    size_type
    FiniteGrid<R>::index(const array<size_type>& pos) const 
    {
      size_type result=0;
      for(dimension_type i=0; i!=dimension(); ++i) {
	result += (_strides[i] * pos[i]);
      }
      return result;
    }

  }
}

#endif /* _GRID_H */
