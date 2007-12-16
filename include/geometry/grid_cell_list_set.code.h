/***************************************************************************
 *            grid_approximation.code.h
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

#include "list_set.h"
#include "grid_cell.h"
#include "grid_block.h"
#include "grid_mask_set.h"
#include "grid_approximation.h"

#include "output/logging.h"




namespace Ariadne {






template<class R>
Geometry::GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
  : _grid_ptr(new Grid<R>(g)), _lattice_set(g.dimension())
{
}


template<class R>
Geometry::GridCellListSet<R>::GridCellListSet(const Grid<R>& g, 
                                    const Combinatoric::LatticeCellListSet& lcls)
  : _grid_ptr(new Grid<R>(g)), _lattice_set(lcls)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(g,lcls,"GridCellListSet::GridCellListSet(Grid g, LatticeCellListSet lcls)");
}


template<class R>
Geometry::GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
  : _grid_ptr(gms._grid_ptr), _lattice_set(gms.dimension())
{
  this->_lattice_set.adjoin(gms._lattice_set);
}


template<class R>
Geometry::GridCellListSet<R>::GridCellListSet(const GridCellListSet<R>& gcls)
  : _grid_ptr(gcls._grid_ptr), _lattice_set(gcls._lattice_set)
{
}

template<class R>
Geometry::GridCellListSet<R>&
Geometry::GridCellListSet<R>::operator=(const GridCellListSet<R>& gcls)
{
  if(this!=&gcls) {
    this->_grid_ptr = gcls._grid_ptr;
    this->_lattice_set=gcls._lattice_set;
  }
  return *this;
}



template<class R>
Geometry::GridCellListSet<R>::operator ListSet< Rectangle<R> >() const
{
  ListSet< Rectangle<R> > result(dimension());
  for(size_type i=0; i!=size(); ++i) {
    result.push_back((*this)[i]);
  }
  return result;
}


template<class R>
Geometry::GridCellListSet<R>*
Geometry::GridCellListSet<R>::clone() const
{
  return new GridCellListSet<R>(*this);
}


template<class R>
tribool
Geometry::GridCellListSet<R>::contains(const Point<R>& pt) const
{
  return !Geometry::disjoint(*this,Rectangle<R>(pt));
}


template<class R>
tribool
Geometry::GridCellListSet<R>::superset(const Rectangle<R>& r) const
{
  return Geometry::subset(r,*this);
}


template<class R>
tribool
Geometry::GridCellListSet<R>::intersects(const Rectangle<R>& r) const
{
  return !Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::GridCellListSet<R>::disjoint(const Rectangle<R>& r) const
{
  return Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::GridCellListSet<R>::subset(const Rectangle<R>& r) const
{
  return Geometry::subset(*this,r);
}


template<class R> 
Geometry::Rectangle<R> 
Geometry::GridCellListSet<R>::bounding_box() const 
{
  return GridBlock<R>(this->grid(),this->lattice_set().bounding_block()); 
}




template<class R>
void
Geometry::GridCellListSet<R>::clear()
{
  this->_lattice_set.clear();
}



template<class R> inline
void 
Geometry::GridCellListSet<R>::restrict_outer_approximation(const SetInterface<R>& s)
{
  Geometry::GridCellListSet<R> result(this->grid());
  Rectangle<R> cell(this->dimension());
  for(typename GridCellListSet<R>::const_iterator cell_iter=this->begin();
      cell_iter!=this->end(); ++cell_iter)
  {
    cell=*cell_iter;
    if(possibly(s.intersects(cell))) {
      result.adjoin(*cell_iter);
    }
  }
  *this=result;
}


template<class R> inline
void 
Geometry::GridCellListSet<R>::restrict_inner_approximation(const SetInterface<R>& s)
{
  Geometry::GridCellListSet<R> result(this->grid());
  Rectangle<R> cell(this->dimension());
  for(typename GridCellListSet<R>::const_iterator cell_iter=this->begin();
      cell_iter!=this->end(); ++cell_iter)
  {
    cell=*cell_iter;
    if(s.superset(cell)) {
      result.adjoin(*cell_iter);
    }
  }
  *this=result;
}





template<class R>
std::ostream&     
Geometry::GridCellListSet<R>::summarize(std::ostream& os) const 
{
  os << "GridCellListSet("
     << " grid=" << this->grid() << ","
     << " size=" << this->size() << ","
     << " )" << std::endl;
  return os;
}

template<class R>
std::ostream&     
Geometry::GridCellListSet<R>::write(std::ostream& os) const 
{
  os << "GridCellListSet("
     << " grid=" << this->grid() << ","
     << " size=" << this->size() << ","
     << " lattice_set=" << this->lattice_set()
     << " )" << std::endl;
  return os;
}


template<class R>
void
Geometry::GridCellListSet<R>::_instantiate()
{
  tribool tb;
  GridBlock<R>* gb=0;
  GridCellListSet<R>* gcls=0;
  
  tb=Geometry::subset(*gcls,*gb);
  tb=Geometry::subset(*gcls,*gcls);
}




// Geometric predicates ---------------------------------------------------

template<class R>
tribool
Geometry::subset(const GridCellListSet<R>& gcls, const GridBlock<R>& gb)
{
  ARIADNE_CHECK_SAME_GRID(gcls,gb,"tribool subset(GridCellListSet gcls, GridBlock gb)");
  return subset(gcls.lattice_set(),gb.lattice_set());
}

template<class R>
tribool
Geometry::subset(const GridCellListSet<R>& gcls1, const GridCellListSet<R>& gcls2)
{
  ARIADNE_CHECK_SAME_GRID(gcls1,gcls2,"tribool subset(GridCellListSet gcls1, GridCellListSet gcls2)");
  return subset(gcls1.lattice_set(),gcls2.lattice_set());
}



} // namespace Ariadne
