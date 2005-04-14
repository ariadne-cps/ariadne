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

#include <boost/numeric/interval.hpp>
#include <binary_word.h>
#include <array.h>

namespace Ariadne {	
namespace Geometry {	
    template<typename R> class Rectangle;

    /*!\brief A grid of rectangles in Euclidean space.
     */
    template<typename R> class Grid {
      private:
	Rectangle<R> _bounding_box;
	array<unsigned int> _subdivision_coordinates;
    };
}
}

#endif /* _GRID_H */
