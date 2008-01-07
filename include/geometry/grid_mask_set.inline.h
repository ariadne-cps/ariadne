/***************************************************************************
 *            grid_mask_set.inline.h
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
 
#include "grid_approximation.h"

namespace Ariadne {


template<class R> inline
bool 
Geometry::GridMaskSet<R>::operator==(const GridMaskSet<R>& gms) 
{
  throw Deprecated(__PRETTY_FUNCTION__);
}


template<class R> inline
bool 
Geometry::GridMaskSet<R>::operator!=(const GridMaskSet<R>& gms) 
{
  return !(*this==gms);
}

template<class R> inline
const Geometry::Grid<R>& 
Geometry::GridMaskSet<R>::grid() const 
{
  return this->_grid; 
}

template<class R> inline
const Combinatoric::LatticeBlock& 
Geometry::GridMaskSet<R>::block() const 
{
  return _lattice_set.block(); 
}

template<class R> inline
const Combinatoric::LatticeMaskSet& 
Geometry::GridMaskSet<R>::lattice_set() const 
{
  return _lattice_set; 
}


template<class R> inline
dimension_type 
Geometry::GridMaskSet<R>::dimension() const 
{
  return _lattice_set.dimension(); 
}

template<class R> inline
size_type 
Geometry::GridMaskSet<R>::capacity() const 
{
  return _lattice_set.capacity(); 
}

//template<class R> inline const SizeArray& GridMaskSet<R>::sizes() const { return _lattice_set.sizes(); } 

template<class R> inline
const BooleanArray& 
Geometry::GridMaskSet<R>::mask() const 
{
  return _lattice_set.mask(); 
}

template<class R> inline
tribool 
Geometry::GridMaskSet<R>::empty() const 
{
  return _lattice_set.empty(); 
}

template<class R> inline
tribool 
Geometry::GridMaskSet<R>::bounded() const 
{
  return _lattice_set.bounded(); 
}

template<class R> inline
size_type 
Geometry::GridMaskSet<R>::size() const 
{
  return _lattice_set.size(); 
}

template<class R> inline
Geometry::GridCell<R> 
Geometry::GridMaskSet<R>::operator[](size_type i) const 
{
  return GridCell<R>(this->_grid,_lattice_set[i]); 
}

template<class R> inline
typename Geometry::GridMaskSet<R>::const_iterator 
Geometry::GridMaskSet<R>::begin() const 
{
  return const_iterator(this->_grid,this->_lattice_set.begin()); 
}
template<class R> inline
typename Geometry::GridMaskSet<R>::const_iterator 
Geometry::GridMaskSet<R>::end() const 
{
  return const_iterator(this->_grid,this->_lattice_set.end()); 
}


template<class R> inline
void 
Geometry::GridMaskSet<R>::remove(const GridCell<R>& gc) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gc,"void GridMaskSet::remove(GridCell gc)");
  this->_lattice_set.remove(gc._lattice_set);
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::remove(const GridCellListSet<R>& gcls) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gcls,"void GridMaskSet::remove(GridCellListSet gcls)");
  this->_lattice_set.remove(gcls._lattice_set);
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::remove(const GridMaskSet<R>& gms) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gms,"void GridMaskSet::remove(GridMaskSet gms)");
  this->_lattice_set.remove(gms._lattice_set);
}

template<class R> inline
void
Geometry::GridMaskSet<R>::adjoin_unbounded_cell() 
{
  this->_lattice_set.adjoin_unbounded_cell();
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::adjoin(const GridCell<R>& gc) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gc,"void GridMaskSet::adjoin(GridCell gc)");
  this->_lattice_set.adjoin(gc._lattice_set);
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::adjoin(const GridBlock<R>& gb) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gb,"void GridMaskSet::adjoin(GridBlock gb)");
  this->_lattice_set.adjoin(gb._lattice_set);
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::adjoin(const GridCellListSet<R>& gcls) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gcls,"void GridMaskSet::adjoin(GridCellListSet gcls)");
  this->_lattice_set.adjoin(gcls._lattice_set);
}


template<class R> inline
void 
Geometry::GridMaskSet<R>::adjoin(const GridMaskSet<R>& gms) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gms,"void GridMaskSet::adjoin(GridMaskSet gms)");
  this->_lattice_set.adjoin(gms._lattice_set);
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::restrict(const GridCellListSet<R>& gcls) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gcls,"void GridMaskSet::restrict(GridCellListSet gcls)");
  this->_lattice_set.restrict(gcls._lattice_set);
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::restrict(const GridMaskSet<R>& gms) 
{
  ARIADNE_CHECK_SAME_GRID(*this,gms,"void GridMaskSet::restrict(GridMaskSet gms)");
  this->_lattice_set.restrict(gms._lattice_set);
}


template<class R> inline
Geometry::GridMaskSet<R> 
Geometry::GridMaskSet<R>::neighbourhood() const 
{
  return GridMaskSet(this->grid(),this->_lattice_set.neighbourhood());
}


template<class R> inline
Geometry::GridMaskSet<R> 
Geometry::GridMaskSet<R>::adjoining() const 
{
  return GridMaskSet(this->grid(),this->_lattice_set.adjoining());
}


template<class R> template<class BS> inline
Geometry::GridMaskSet<R>::operator ListSet<BS> () const 
{
  ListSet<BS> result(this->dimension());
  Box<R> r(this->dimension());
  BS bs(this->dimension());
  for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
    r=*iter;
    bs=BS(r);
    result.push_back(bs);
  }
  return result;
}

template<class R> template<class BS>
void 
Geometry::GridMaskSet<R>::adjoin_outer_approximation(const ListSet<BS>& ls)
{
  for(size_type i=0; i!=ls.size(); ++i) {
    this->adjoin_outer_approximation(ls[i]);
  }
}


template<class R> inline
void 
Geometry::GridMaskSet<R>::adjoin_over_approximation(const Box<R>& r)
{
  this->adjoin(over_approximation(r,this->grid()));
}

template<class R> inline
void 
Geometry::GridMaskSet<R>::adjoin_under_approximation(const Box<R>& r)
{
  this->adjoin(under_approximation(r,this->grid()));
}


template<class R> template<class BS> inline
void 
Geometry::GridMaskSet<R>::adjoin_outer_approximation(const BS& bs)
{
  this->adjoin(outer_approximation(bs,this->grid()));
}


template<class R> template<class BS> inline
void 
Geometry::GridMaskSet<R>::adjoin_inner_approximation(const BS& bs)
{
  this->adjoin(inner_approximation(bs,this->grid()));
}












template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const GridMaskSet<R>& gms) {
  return gms.write(os);
}



}
