/***************************************************************************
 *            grid_cell.code.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <ostream>

#include "../geometry/rectangle.h"
#include "../geometry/grid_cell.h"
#include "../geometry/grid_block.h"

namespace Ariadne {


template<class R> inline
Geometry::GridBlock<R>
Geometry::GridCell<R>::neighbourhood() const 
{
  return GridBlock<R>(this->_grid_ptr,this->_lattice_set.neighbourhood());
}


template<class R>
std::ostream&
Geometry::GridCell<R>::write(std::ostream& os) const 
{
  return os << Rectangle<R>(*this);
}


}
