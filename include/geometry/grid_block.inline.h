/***************************************************************************
 *            grid_block.h
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
 
namespace Ariadne {



template<class R> inline
Geometry::GridBlock<R>::GridBlock(const Grid<R>* gptr, const Combinatoric::LatticeBlock& lc)
  : _grid_ptr(gptr), _lattice_set(lc)
{
}


template<class R> inline
const Geometry::Grid<R>& 
Geometry::GridBlock<R>::grid() const 
{
  return *this->_grid_ptr; 
}


template<class R> inline
dimension_type 
Geometry::GridBlock<R>::dimension() const 
{
  return this->_lattice_set.dimension(); 
}


template<class R> inline
const Combinatoric::LatticeBlock& 
Geometry::GridBlock<R>::lattice_set() const 
{
  return this->_lattice_set; 
}


template<class R> inline
tribool 
Geometry::GridBlock<R>::empty() const 
{
  return this->_lattice_set.empty(); 
}


template<class R> inline
tribool 
Geometry::GridBlock<R>::bounded() const 
{
  return true; 
}


template<class R> inline
Geometry::Box<R> 
Geometry::GridBlock<R>::bounding_box() const 
{
  return *this; 
}



template<class R> inline
typename Geometry::GridBlock<R>::const_iterator 
Geometry::GridBlock<R>::begin() const 
{
  return const_iterator(*this->_grid_ptr,_lattice_set.begin()); 
}


template<class R> inline
typename Geometry::GridBlock<R>::const_iterator 
Geometry::GridBlock<R>::end() const 
{
  return const_iterator(*this->_grid_ptr,_lattice_set.end()); 
}



template<class R> inline
std::ostream&
Geometry::operator<<(std::ostream& os, const GridBlock<R>& gb) {
  return gb.write(os);
}


}
