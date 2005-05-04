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

namespace Ariadne {	
namespace Geometry {	
    template<typename R> class Rectangle;
    template<typename R> class State;
    
    /*! \brief A rectangle defined on a grid. 
     */
    template<typename R>
    class GridRectangle {
      public:
	typedef R real_type;
	typedef S state_type;
	typedef size_t size_type;
	
	/*!\brief Construct from a grid and a binary word.
	 */
	GridRectangle(const Grid& g, const BinaryWord& w) 
	    : _bounding_box(g.bounding_box()), 
	      _subdivision_coordinates(g.subdivision_coordinates()), 
	      _word(w) 
	{ } 

	/*!\brief Convert to an ordinary rectangle.
	 */
	operator Rectangle<R>() const; 	    
      private:
	const Rectangle& _bounding_box;
	const array<size_type> & _subdivision_coordinates;
	BinaryWord _word;
    };

    inline GridRectangle::operator Rectangle<R>() const {
	Rectangle res(_bounding_box);
	array<size_type>::const_iterator coord_iter=_subdivision_coordinates.begin();
	BinaryWord::const_iterator word_iter=_word.begin();
	
	while(word_iter!=_word.end()) {
	    size_type i=(*coord_iter);
	    R centre = ( res.lower(i) + res.upper(i) ) / R(2);
	    if( (*word_iter)==0 ) {
		res.upper(i)=centre; }
	    else {
		res.lower(i)=centre;
	    }
	}
	
	return res;
    }
		
}
}

#endif /* _GRID_RECTANGLE_H */
