/***************************************************************************
 *            grid_rectangle.h
 *
 *  18 January 2005
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
 
/*! \file grid_rectangle.h
 *  \brief Rectangles defined on a grid.
 */

#ifndef _GRID_RECTANGLE_H
#define _GRID_RECTANGLE_H

#include "interval.h"
#include "binary_word.h"
#include "rectangle.h"
#include "state.h"
#include "grid.h"

namespace Ariadne {	
  namespace Geometry {	
    template<typename R> class Rectangle;
    template<typename R> class State;
    
    /*! \brief A rectangle defined on a grid. 
     */
    template<typename R>
    class FiniteGridCell {
     public:
      typedef R real_type;
      typedef State<R> state_type;
   
      /*!\brief Construct from a grid and an integer array.
       */
      FiniteGridCell(const FiniteGrid<R>& g, 
		       const array<size_type>& pos) 
	: _grid(g), _position(pos) { } 

      dimension_type dimension() const { return _grid.dimension(); }

      /*!\brief Convert to an ordinary rectangle.
       */
      operator Rectangle<R>() const; 	    
     private:
      const FiniteGrid<R>& _grid;
      array<size_type> _position;
    };


    /*! \brief A rectangle defined on a grid. 
     */
    template<typename R>
    class FiniteGridRectangle {
     public:
      typedef R real_type;
      typedef State<R> state_type;
   
      /*!\brief Construct from a grid and an integer array.
       */
      FiniteGridRectangle(const FiniteGrid<R>& g, 
			  const array<size_type>& l,
			  const array<size_type>& u)
	: _grid(g), _lower(l), _upper(u) 
      {
	assert(_lower.size()==_upper.size());
	for(size_type i=0; i!=_lower.size(); ++i) {
	  assert(_lower[i]!=_upper[i]);
	  if(_lower[i]>_upper[i]) {
	    std::swap(_lower[i],_upper[i]);
	  }
	  assert(_upper[i] < g.size(i));
	}
      }

      dimension_type dimension() const { return _grid.dimension(); }

      /*!\brief Convert to an ordinary rectangle.
       */
      operator Rectangle<R>() const; 	    
     private:
      const FiniteGrid<R>& _grid;
      array<size_type> _lower;
      array<size_type> _upper; 
   };

    /*! \brief A rectangle defined on a grid. 
     */
    template<typename R>
    class InfiniteGridCell {
     public:
      typedef R real_type;
      typedef State<R> state_type;
	
      /*!\brief Construct from a grid and an integer array.
       */
      InfiniteGridCell(const InfiniteGrid<R>& g, 
		       const array<int>& pos) 
	: _grid(g), _position(pos) { } 

      dimension_type dimension() const { return _grid.dimension(); }

      /*!\brief Convert to an ordinary rectangle.
       */
      operator Rectangle<R>() const; 	    
     private:
      const InfiniteGrid<R>& _grid;
      array<int> _position;
    };



    /*! \brief A rectangle defined on a grid. 
     */
    template<typename R>
    class InfiniteGridRectangle {
     public:
      typedef R real_type;
      typedef State<R> state_type;
	
      /*!\brief Construct from a grid and an integer array.
       */
      InfiniteGridRectangle(const InfiniteGrid<R>& g, 
		       const array<int>& l,
		       const array<int>& u) 
	: _grid(g), _lower(l), _upper(u) 
      { 
	assert(_lower.size() == _upper.size());
	for(size_type i=0; i!=_lower.size(); ++i) {
	  assert(_lower[i] != _upper[i]);
	  if(_lower[i] > _upper[i]) {
	    std::swap(_lower[i],_upper[i]);
	  }
	}
      } 
      
      dimension_type dimension() const { return _grid.dimension(); }

      /*!\brief Convert to an ordinary rectangle.
       */
      operator Rectangle<R>() const; 	    
     private:
      const InfiniteGrid<R>& _grid;
      array<int> _lower;
      array<int> _upper;
    };



    /*! \brief A rectangle defined on a grid. 
     */
    template<typename R>
    class PartitionGridCell {
     public:
      typedef R real_type;
      typedef State<R> state_type;
      typedef size_t size_type;
	
      /*!\brief Construct from a grid and a binary word.
       */
      PartitionGridCell(const PartitionGrid<R>& g, 
			const BinaryWord& w) 
	: _grid(g), _word(w) 
      { } 

      dimension_type dimension() const { return _grid.dimension(); }

      /*!\brief Convert to an ordinary rectangle.
       */
      operator Rectangle<R>() const; 	    
     private:
      const PartitionGrid<R>& _grid;
      BinaryWord _word;
    };


    template<typename R>
    inline 
    InfiniteGridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());
      
      for(dimension_type i=0; i!=dimension(); ++i) {
	result.set_lower(i, _grid.subdivision_length(i) * _position[i]);
	result.set_upper(i, _grid.subdivision_length(i) * (_position[i] + 1));
      }
      
      return result;
    }


    template<typename R>
    inline 
    InfiniteGridRectangle<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());
      
      for(size_type i=0; i!=dimension(); ++i) {
	result.set_lower(i, _grid.subdivision_position(i,_lower[i]));
	result.set_upper(i, _grid.subdivision_position(i,_upper[i]));
      }
      
      return result;
    }


    template<typename R>
    inline 
    FiniteGridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(_grid.dimension());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
	result.set_lower(i,_grid.subdivision_position(i,_position[i]));
	result.set_upper(i,_grid.subdivision_position(i,_position[i]+1));
      }

      return result;
    }

    
    template<typename R>
    inline 
    FiniteGridRectangle<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
	result.set_lower(i, _grid.subdivision_position(i,_lower[i]));
	result.set_upper(i, _grid.subdivision_position(i,_upper[i]));
      }

      return result;
    }

    
    template<typename R>
    inline 
    PartitionGridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> res(_grid.bounding_box());
      sequence<unsigned short>::const_iterator coord_iter=_grid.subdivision_coordinates().begin();
      BinaryWord::const_iterator word_iter=_word.begin();
      
      while(word_iter!=_word.end()) {
	size_type i=(*coord_iter);
	R centre = ( res.lower(i) + res.upper(i) ) / R(2);
	if( (*word_iter)==0 ) {
	  res.set_upper(i,centre); }
	else {
	  res.set_lower(i,centre);
	}
      }
      
      return res;
    }

    
  }
}

#endif /* _GRID_RECTANGLE_H */
