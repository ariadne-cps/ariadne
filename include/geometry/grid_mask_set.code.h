/***************************************************************************
 *            grid_mask_set.code.h
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

#include "grid_mask_set.h"

#include <ostream>
#include <sstream>
#include <string>

#include "base/stlio.h"

#include "combinatoric/array_operations.h"

#include "geometry/box.h"

#include "geometry/grid_cell.h"
#include "geometry/grid_block.h"
#include "geometry/grid_cell_list_set.h"

#include "geometry/list_set.h"
#include "geometry/partition_tree_set.h"

#include "geometry/set_interface.h"

#include "output/logging.h"





namespace Ariadne {


template<class R>
GridMaskSet<R>::GridMaskSet()
  : _grid(), _lattice_set(0) 
{ 
}

template<class R>
GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg)
  : _grid(fg.grid()), _lattice_set(fg.lattice_block()) 
{ 
}


template<class R>
GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg, const BooleanArray& m)
  : _grid(fg.grid()), _lattice_set(fg.lattice_block(),m)
{
}


template<class R>
GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Box<R>& bb)
  : _grid(g), _lattice_set(g.index_block(bb)) 
{ 
}

template<class R>
GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeBlock& b)
  : _grid(g), _lattice_set(b) 
{ 
}


template<class R>
GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeBlock& b, const BooleanArray& m)
  : _grid(g), _lattice_set(b,m)
{
}


template<class R>
GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeMaskSet& ms)
  : _grid(g), _lattice_set(ms)
{
}


template<class R>
GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
  : _grid(gms._grid), _lattice_set(gms._lattice_set)
{
}


template<class R>
GridMaskSet<R>&
GridMaskSet<R>::operator=(const GridMaskSet<R>& gms) 
{
  if(this!=&gms) {
    this->_grid=gms._grid;
    this->_lattice_set=gms._lattice_set;
  }
  return *this;
}


template<class R>
GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
  : _grid(gcls._grid),
    _lattice_set(gcls.lattice_set())
{
}



template<class R>
FiniteGrid<R>
GridMaskSet<R>::finite_grid() const
{
  return FiniteGrid<R>(this->grid(),this->block());
}


template<class R>
Box<R>
GridMaskSet<R>::extent() const
{
  return this->bounds();
}


template<class R> inline
GridBlock<R> 
GridMaskSet<R>::bounds() const 
{
  return GridBlock<R>(this->_grid,_lattice_set.block()); 
}


// FIXME: Memory leak

template<class R>
GridMaskSet<R>*
GridMaskSet<R>::clone() const
{
  return new GridMaskSet<R>(*this);
}


template<class R>
tribool
GridMaskSet<R>::contains(const Point<R>& pt) const
{
  return !Ariadne::disjoint(*this,Box<R>(pt));
}


template<class R>
tribool
GridMaskSet<R>::superset(const Box<R>& r) const
{
  return Ariadne::subset(r,*this);
}


template<class R>
tribool
GridMaskSet<R>::intersects(const Box<R>& r) const
{
  return !Ariadne::disjoint(*this,r);
}


template<class R>
tribool
GridMaskSet<R>::disjoint(const Box<R>& r) const
{
  return Ariadne::disjoint(*this,r);
}


template<class R>
tribool
GridMaskSet<R>::subset(const Box<R>& r) const
{
  return Ariadne::subset(*this,r);
}


template<class R> 
Box<R> 
GridMaskSet<R>::bounding_box() const 
{
  return GridBlock<R>(grid(),bounds()); 
}


template<class R>
void
GridMaskSet<R>::clear()
{
  this->_lattice_set.clear();
}




template<class R>
void
GridMaskSet<R>::_instantiate_geometry_operators()
{
  typedef Interval<R> I;
  tribool tb;
  Box<R>* r=0;
  
  //FiniteGrid<R>* fg=0;
  GridCell<R>* gc=0;
  GridBlock<R>* gb=0;
  GridCellListSet<R>* gcls=0;
  GridMaskSet<R>* gms=0;
  
  tb=Ariadne::subset(*r,*gms);
  tb=Ariadne::subset(*gms,*r);
  tb=Ariadne::superset(*gms,*r);
  tb=Ariadne::disjoint(*r,*gms);
  tb=Ariadne::disjoint(*gms,*r);
  
  tb=Ariadne::overlap(*gb,*gms);
  tb=Ariadne::overlap(*gms,*gb);
  tb=Ariadne::overlap(*gcls,*gms);
  tb=Ariadne::overlap(*gms,*gcls);
  tb=Ariadne::overlap(*gms,*gms);
  
  tb=Ariadne::subset(*gc,*gms);
  tb=Ariadne::subset(*gb,*gms);
  tb=Ariadne::subset(*gcls,*gms);
  tb=Ariadne::subset(*gms,*gms);
  
  *gms=Ariadne::regular_intersection(*gb,*gms);
  *gms=Ariadne::regular_intersection(*gms,*gb);
  *gcls=Ariadne::regular_intersection(*gcls,*gms);
  *gcls=Ariadne::regular_intersection(*gms,*gcls);
  *gms=Ariadne::regular_intersection(*gms,*gms);
  *gcls=Ariadne::difference(*gcls,*gms);
  *gms=Ariadne::difference(*gms,*gms);
  *gms=Ariadne::join(*gms,*gms);
}




// Geometric predicates ---------------------------------------------------



template<class R>
tribool
disjoint(const GridBlock<R>& gb, const GridMaskSet<R>& gms) {
  ARIADNE_CHECK_SAME_GRID(gb,gms,"tribool disjoint(GridBlock gb, GridMaskSet gms)");
  return disjoint(gb.lattice_set(),gms.lattice_set());
}


template<class R>
tribool
disjoint(const GridMaskSet<R>& gms, const GridBlock<R>& gb) {
  ARIADNE_CHECK_SAME_GRID(gms,gb,"tribool disjoint(GridMaskSet gms, GridBlock gb)");
  return disjoint(gms.lattice_set(),gb.lattice_set());
}


template<class R>
tribool
disjoint(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"tribool disjoint(GridMaskSet gms1, GridMaskSet gms2)");
  return disjoint(gms1.lattice_set(),gms2.lattice_set());
}


template<class R>
tribool
disjoint(const Box<R>& r, const GridMaskSet<R>& gms) 
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gms,"tribool disjoint(Box r, GridMaskSet gms)");
  Box<R> br=closed_intersection(r,Box<R>(gms.bounding_box()));
  GridBlock<R> gb=outer_approximation(br,gms.grid());
  return !overlap(gb,gms);
}


template<class R>
tribool
disjoint(const GridMaskSet<R>& gms, const Box<R>& r) {
  return disjoint(r,gms);
}





template<class R>
tribool
overlap(const GridBlock<R>& gb, const GridMaskSet<R>& gms) {
  ARIADNE_CHECK_SAME_GRID(gb,gms,"tribool overlap(GridBlock gb, GridMaskSet gms)");
  return overlap(gb.lattice_set(),gms.lattice_set());
}


template<class R>
tribool
overlap(const GridMaskSet<R>& gms, const GridBlock<R>& gb) {
  ARIADNE_CHECK_SAME_GRID(gms,gb,"tribool overlap(GridMaskSet gms, GridBlock gb)");
  return overlap(gms.lattice_set(),gb.lattice_set());
}


template<class R>
tribool
overlap(const GridCellListSet<R>& A, const GridMaskSet<R>& B) {
  ARIADNE_CHECK_SAME_GRID(A,B,"overlap(GridCellListSet<R>,GridMaskSet<R>)");
  return overlap(A.lattice_set(),B.lattice_set());
}


template<class R>
tribool
overlap(const GridMaskSet<R>& A, const GridCellListSet<R>& B) {
  ARIADNE_CHECK_SAME_GRID(A,B,"overlap(GridMaskSet<R>,GridCellListSet<R>)");
  return overlap(A.lattice_set(),B.lattice_set());
}


template<class R>
tribool
overlap(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"tribool overlap(GridMaskSet gms2, GridMaskSet gms2)");
  return overlap(gms1.lattice_set(),gms2.lattice_set());
}






template<class R>
tribool
subset(const GridMaskSet<R>& gms, const GridBlock<R>& gb)
{
  ARIADNE_CHECK_SAME_GRID(gms,gb,"tribool subset(GridMaskSet gms, GridBlock gb)");
  return subset(gms.lattice_set(),gb.lattice_set());
}

template<class R>
tribool
subset(const GridCell<R>& gc, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gc,gms,"tribool subset(GridCell gc, GridMaskSet gms)");
  return subset(gc.lattice_set(),gms.lattice_set()); 
}

template<class R>
tribool
subset(const GridBlock<R>& gb, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gb,gms,"tribool subset(GridBlock gb, GridMaskSet gms)");
  return subset(gb.lattice_set(),gms.lattice_set()); 
}

template<class R>
tribool
subset(const GridCellListSet<R>& gcls, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gcls,gms,"tribool subset(GridCellListSet gcls, GridMaskSet gms)");
  return subset(gcls.lattice_set(),gms.lattice_set()); 
}

template<class R>
tribool
subset(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"tribool subset(GridMaskSet gms1, GridMaskSet gms2)");
  return subset(gms1.lattice_set(),gms2.lattice_set());
}


template<class R>
tribool
superset(const GridMaskSet<R>& gms, const Box<R>& r)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gms,"tribool superset(GridMaskSet gms, Box r)");
  if(!subset(r,gms.bounding_box())) {
    return false;
  }
  return subset(outer_approximation(r,gms.grid()),gms);
}



template<class R>
tribool
subset(const Box<R>& r, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gms,"tribool subset(Box r, GridMaskSet gms)");
  if(!subset(r,gms.bounding_box())) {
    return false;
  }
  return subset(outer_approximation(r,gms.grid()),gms);
}


template<class R>
tribool
subset(const GridMaskSet<R>& gms, const Box<R>& r)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(gms,r,"tribool subset(GridMaskSet gms, Box r)");
  if(!subset(gms,r.bounding_box())) {
    return false;
  }
  return subset(gms,outer_approximation(r,gms.grid()));
}





template<class R>
GridMaskSet<R>
join(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"GridMaskSet join(GridMaskSet gms1, GridMaskSet gms2)");
  if(gms1.block()==gms2.block()) { 
    return GridMaskSet<R>(gms1.grid(), join(gms1.lattice_set(),gms2.lattice_set()));
  } else {
    GridMaskSet<R> result(gms1.grid(),rectangular_hull(gms1.block(),gms2.block()));
    result.adjoin(gms1);
    result.adjoin(gms2);
    return result;
  }
}


template<class R>
GridMaskSet<R>
regular_intersection(const GridMaskSet<R>& gms, const GridBlock<R>& gb)
{
  ARIADNE_CHECK_SAME_GRID(gms,gb,"GridMaskSet regular_intersection(GridMaskSet gms, GridBlock gb)");
  return GridMaskSet<R>(gms.grid(), regular_intersection(gms.lattice_set(),gb.lattice_set()));
}

template<class R>
GridMaskSet<R>
regular_intersection(const GridBlock<R>& gb, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gb,gms,"GridMaskSet regular_intersection(GridBlock gb, GridMaskSet gms)");
  return GridMaskSet<R>(gb.grid(), regular_intersection(gb.lattice_set(),gms.lattice_set()));
}

template<class R>
GridMaskSet<R>
regular_intersection(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"GridMaskSet regular_intersection(GridMaskSet gms1, GridMaskSet gms2)");
  if(gms1.block()==gms2.block()) { 
    return GridMaskSet<R>(gms1.grid(), regular_intersection(gms1.lattice_set(),gms2.lattice_set()));
  } else {
    GridMaskSet<R> result(gms1.grid(),rectangular_hull(gms1.block(),gms2.block()));
    result.adjoin(gms1);
    result.restrict(gms2);
    return result;
  }
}

template<class R>
GridCellListSet<R>
regular_intersection(const GridCellListSet<R>& gcls, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gcls,gms,"GridCellListSet regular_intersection(GridCellListSet gcls, GridMaskSet gms)");
  return GridCellListSet<R>(gcls.grid(), regular_intersection(gcls.lattice_set(),gms.lattice_set()));
}

template<class R>
GridCellListSet<R>
regular_intersection(const GridMaskSet<R>& gms, const GridCellListSet<R>& B)
{
  ARIADNE_CHECK_SAME_GRID(gms,B,"GridCellListSet regular_intersection(GridMaskSet<R>,GridCellListSet<R>)");
  return GridCellListSet<R>(gms.grid(), regular_intersection(gms.lattice_set(),B.lattice_set()));
}


template<class R>
GridMaskSet<R>
difference(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"GridMaskSet difference(GridMaskSet<R>,GridMaskSet<R>)");
  if(gms1.block()==gms2.block()) { 
    return GridMaskSet<R>(gms1.grid(), difference(gms1.lattice_set(),gms2.lattice_set()));
  } else {
    GridMaskSet<R> result(gms1.grid(),rectangular_hull(gms1.block(),gms2.block()));
    result.adjoin(gms1);
    result.remove(gms2);
    return result;
  }
}

template<class R>
GridCellListSet<R>
difference(const GridCellListSet<R>& gcls, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gcls,gms,"difference(GridCellListSet gcls, GridMaskSet gms)");
  return GridCellListSet<R>(gcls.grid(),difference(gcls.lattice_set(),gms.lattice_set()));
}





template<class R> 
void 
GridMaskSet<R>::adjoin_over_approximation(const Box<R>& r)
{
  this->adjoin(over_approximation(r,this->grid()));
}

template<class R> 
void 
GridMaskSet<R>::adjoin_under_approximation(const Box<R>& r)
{
  this->adjoin(under_approximation(r,this->grid()));
}

template<class R> 
void 
GridMaskSet<R>::adjoin_outer_approximation(const SetInterface< Box<R> >& s)
{
  FiniteGrid<R> fg(this->grid(),this->block());
  this->adjoin(outer_approximation(s,fg));
}

template<class R> 
void 
GridMaskSet<R>::adjoin_inner_approximation(const SetInterface< Box<R> >& s)
{
  FiniteGrid<R> fg(this->grid(),this->block());
  this->adjoin(inner_approximation(s,fg));
}

template<class R> 
void 
GridMaskSet<R>::restrict_outer_approximation(const SetInterface< Box<R> >& s)
{
  FiniteGrid<R> fg(this->grid(),this->block());
  this->restrict(outer_approximation(s,fg));
}

template<class R> 
void 
GridMaskSet<R>::restrict_inner_approximation(const SetInterface< Box<R> >& s)
{
  FiniteGrid<R> fg(this->grid(),this->block());
  this->restrict(inner_approximation(s,fg));
}






template<class R>
std::string 
GridMaskSet<R>::summary() const 
{
  std::stringstream ss;
  ss << "GridMaskSet( "
     << " grid=" << this->grid() << ","
     << " block=" << this->block() << ","
     << " extent=" << Box<R>(this->bounds()) << ","
     << " size=" << this->size() << ","
     << " capacity=" << this->capacity() << ","
     << " )";
  return ss.str();
}

template<class R>
std::ostream& 
GridMaskSet<R>::write(std::ostream& os) const 
{
  os << "GridMaskSet( " << std::flush;
  os << " grid=" << this->grid() << ",";
  os << " block=" << this->block() << ",";
  os << " extent=" << Box<R>(this->bounds()) << ",";
  os << " size=" << this->size() << ",";
  os << " capacity=" << this->capacity() << ",";
  os << " mask=" << this->mask() << std::endl;
  os << " )";
  return os;
}


} // namespace Ariadne
                                                    
