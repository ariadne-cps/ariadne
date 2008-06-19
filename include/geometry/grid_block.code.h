/***************************************************************************
 *            grid_block.code.h
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

#include "grid_block.h"

#include <ostream>

#include "base/stlio.h"

#include "combinatoric/array_operations.h"

#include "geometry/box.h"

#include "output/logging.h"




namespace Ariadne {
    
template<class R>
GridBlock<R>::GridBlock(const Grid<R>& g)
  : _grid(g), _lattice_set(g.dimension())
{ 
  _lattice_set.set_lower_bound(0,1);
  //_lattice_set.set_lower_bound(0,0);
  _lattice_set.set_upper_bound(0,0);
}




template<class R>
GridBlock<R>::GridBlock(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
  : _grid(g), _lattice_set(l,u)
{
  ARIADNE_CHECK_DIMENSION(g,l.size(),"GridBlock::GridBlock(Grid g, IndexArray l, IndexArray u)");
}


template<class R>
GridBlock<R>::GridBlock(const Grid<R>& g, const Box<R>& r)
  : _grid(g), _lattice_set(g.dimension())
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(g,r,"GridBlock::GridBlock(Grid g,Box r)");
  for(dimension_type i=0; i!=dimension(); ++i) {
    /* TODO: Catch and rethrow exceptions */
    _lattice_set.set_lower_bound(i,g.subdivision_index(i,r.lower_bound(i)));
    _lattice_set.set_upper_bound(i,g.subdivision_index(i,r.upper_bound(i)));
  }
}


template<class R>
GridBlock<R>::GridBlock(const GridCell<R>& gc)
  : _grid(gc._grid), _lattice_set(gc.lattice_set())
{
}


template<class R>
GridBlock<R>::GridBlock(const GridBlock<R>& gb)
  : _grid(gb._grid), _lattice_set(gb.lattice_set())
{
}


template<class R>
GridBlock<R>&
GridBlock<R>::operator=(const GridBlock<R>& gb)
{
  if(this!=&gb) {
    this->_grid = gb._grid;
    this->_lattice_set=gb._lattice_set;
  }
  return *this;
}


template<class R>
R
GridBlock<R>::lower_bound(dimension_type i) const 
{
  return _grid.subdivision_coordinate(i,_lattice_set.lower_bound(i));
}


template<class R>
R
GridBlock<R>::upper_bound(dimension_type i) const 
{
  return _grid.subdivision_coordinate(i,_lattice_set.upper_bound(i));
}

template<class R>
GridBlock<R>
GridBlock<R>::neighbourhood() const 
{
  return GridBlock<R>(this->_grid,this->_lattice_set.neighbourhood());
}


template<class R>
void
GridBlock<R>::_instantiate_geometry_operators() 
{
  typedef Interval<R> I;
  tribool tb;
  Box<R>* r=0;
  //Grid<R>* g=0;
  GridCell<R>* gc=0;
  GridBlock<R>* gb=0;
  
  tb=subset(*r,*gb);
  
  tb=overlap(*gb,*gb);
  tb=subset(*gc,*gb);
  tb=subset(*gb,*gb);
}





// Geometric predicates ---------------------------------------------------

template<class R>
tribool
disjoint(const GridBlock<R>& gb1, const GridBlock<R>& gb2) {
  ARIADNE_CHECK_SAME_GRID(gb1,gb2,"tribool disjoint(GridBlock bg1, GridBlock gb2)");
  return disjoint(gb1.lattice_set(),gb2.lattice_set());
}





template<class R>
tribool
overlap(const GridBlock<R>& gb1, const GridBlock<R>& gb2) {
  ARIADNE_CHECK_SAME_GRID(gb1,gb2,"tribool overlap(GridBlock gb1,GridBlock gb2)");
  return overlap(gb1.lattice_set(),gb2.lattice_set());
}






template<class R>
tribool
subset(const GridCell<R>& gc, const GridBlock<R>& gb)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(gc,gb,"tribool subset(GridCell gc, GridBlock gb)");
  if(gc.grid()==gb.grid()) {
    return subset(gc.lattice_set(),gb.lattice_set());
  }
  return subset(Box<R>(gc),Box<R>(gb));
}

template<class R>
tribool
subset(const GridBlock<R>& gb1, const GridBlock<R>& gb2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(gb1,gb2,"tribool subset(GridBlock gb1, GridBlock gb2)");
  if(gb1.grid()==gb2.grid()) {
    return subset(gb1.lattice_set(),gb2.lattice_set());
  }
  return subset(Box<R>(gb1),Box<R>(gb2));
}


template<class R>
tribool
subset(const Box<R>& r, const GridBlock<R>& gb)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gb,"tribool subset(Box r, GridBlock gb)");
  return subset(r,Box<R>(gb));
}












// Input/output------------------------------------------------------------


template<class R>
std::ostream&
GridBlock<R>::write(std::ostream& os) const 
{
  os << "GridBlock(\n";
  os << "  grid=" << this->grid() << ",\n";
  os << "  size=" << this->lattice_set().size() << ",\n";
  os << "  lattice_set=" << this->lattice_set();
  os << "  rectangle=" << Box<R>(*this);
  os << ")" << std::endl;
  return os;
}




}
