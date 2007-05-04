/***************************************************************************
 *            grid_cell_list_set.inline.h
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
const Geometry::Grid<R>& 
Geometry::GridCellListSet<R>::grid() const 
{
  return *this->_grid_ptr; 
}

template<class R> inline
dimension_type 
Geometry::GridCellListSet<R>::dimension() const 
{
  return this->_lattice_set.dimension(); 
}

template<class R> inline
tribool 
Geometry::GridCellListSet<R>::empty() const 
{
  return this->_lattice_set.empty(); 
}

template<class R> inline
tribool 
Geometry::GridCellListSet<R>::bounded() const 
{
  return true; 
}

template<class R> inline
size_type 
Geometry::GridCellListSet<R>::size() const 
{
  return _lattice_set.size(); 
}

template<class R> inline
const Combinatoric::LatticeCellListSet& 
Geometry::GridCellListSet<R>::lattice_set() const 
{
  return _lattice_set; 
}

template<class R> inline
Geometry::GridCell<R> 
Geometry::GridCellListSet<R>::operator[] (const size_type i) const 
{
  return GridCell<R>(grid(),_lattice_set[i]); 
}


template<class R> inline
typename Geometry::GridCellListSet<R>::const_iterator 
Geometry::GridCellListSet<R>::begin() const 
{
  return const_iterator(*this->_grid_ptr,_lattice_set.begin()); 
}


template<class R> inline
typename Geometry::GridCellListSet<R>::const_iterator 
Geometry::GridCellListSet<R>::end() const 
{
  return const_iterator(*this->_grid_ptr,_lattice_set.end()); 
}


template<class R> inline
void
Geometry::GridCellListSet<R>::unique_sort()
{
  this->_lattice_set.unique_sort();
}


template<class R> inline
void 
Geometry::GridCellListSet<R>::adjoin(const GridCell<R>& c) 
{
  _lattice_set.adjoin(c.lattice_set()); 
}


template<class R> inline
void 
Geometry::GridCellListSet<R>::adjoin(const GridBlock<R>& bl) 
{
  _lattice_set.adjoin(bl.lattice_set()); 
}


template<class R> inline
void 
Geometry::GridCellListSet<R>::adjoin(const GridCellListSet<R>& cls) 
{
  _lattice_set.adjoin(cls.lattice_set()); 
}


template<class R> inline
void 
Geometry::GridCellListSet<R>::adjoin_over_approximation(const Rectangle<R>& r) 
{
  this->adjoin(over_approximation(r,this->grid()));
}


template<class R> template<class BS> inline
void 
Geometry::GridCellListSet<R>::adjoin_outer_approximation(const BS& bs)
{
  this->adjoin(outer_approximation(bs,this->grid()));
}

template<class R> template<class BS> inline
void 
Geometry::GridCellListSet<R>::adjoin_inner_approximation(const BS& bs)
{
  this->adjoin(inner_approximation(bs,this->grid()));
}








template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const GridCellListSet<R>& gcls) {
  return gcls.write(os);
}




}
