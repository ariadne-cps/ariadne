/***************************************************************************
 *            grid_cell.inline.h
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
 
#include "geometry/rectangle.h"

namespace Ariadne {

template<class R> inline 
const Geometry::Grid<R>& 
Geometry::GridCell<R>::grid() const 
{ 
  return *this->_grid_ptr; 
}

template<class R> inline
dimension_type 
Geometry::GridCell<R>::dimension() const 
{
  return this->_lattice_set.dimension(); 
}

template<class R> inline
const Combinatoric::LatticeCell& 
Geometry::GridCell<R>::lattice_set() const 
{
  return this->_lattice_set; 
}

template<class R> inline
tribool 
Geometry::GridCell<R>::bounded() const { 
  return true; 
}

template<class R> inline
Geometry::Rectangle<R> 
Geometry::GridCell<R>::bounding_box() const 
{
  return *this; 
}


template<class R> inline
Geometry::GridCell<R>::GridCell(const Grid<R>& g, const Combinatoric::LatticeCell& pos)
  : _grid_ptr(&g), _lattice_set(pos)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(g,pos,"GridCell<R>::GridCell(Grid<R>,LatticeCell");
}


template<class R> inline
Geometry::GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
  : _grid_ptr(&g), _lattice_set(pos)
{
  ARIADNE_CHECK_DIMENSION(g,pos.size(),"GridCell::GridCell(Grid,IndexArray)");
}


template<class R> inline
R
Geometry::GridCell<R>::lower_bound(dimension_type i) const 
{
  return _grid_ptr->subdivision_coordinate(i,_lattice_set.lower_bound(i));
}


template<class R> inline
R
Geometry::GridCell<R>::upper_bound(dimension_type i) const 
{
  return _grid_ptr->subdivision_coordinate(i,_lattice_set.upper_bound(i));
}


template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const GridCell<R>& gc) {
  return gc.write(os);
}
    

}
