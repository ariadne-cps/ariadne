/***************************************************************************
 *            grid_cell_list_set.code.h
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

#include "grid_cell_list_set.h"

#include <ostream>

#include "../base/stlio.h"

#include "../combinatoric/array_operations.h"

#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/grid_cell.h"
#include "../geometry/grid_block.h"
#include "../geometry/grid_mask_set.h"
#include "../geometry/list_set.h"
#include "../geometry/partition_tree_set.h"

#include "../geometry/set_interface.h"
#include "../geometry/set_reference.h"

#include "../output/logging.h"



namespace {

using namespace Ariadne;


template<class R, class BS>
inline
Geometry::GridCellListSet<R>
outer_approximation_of_basic_set(const BS& bs, const Geometry::Grid<R>& g) 
{
  Geometry::GridCellListSet<R> gcls(g);
  Geometry::GridBlock<R> gbb=outer_approximation(bs.bounding_box(),g);
  Geometry::Rectangle<R> r(bs.dimension());
  //std::cerr << "bs="<<bs<<" bb="<<bs.bounding_box()<<" gbb="<<gbb<<std::endl;
  for(typename Geometry::GridBlock<R>::const_iterator iter=gbb.begin(); iter!=gbb.end(); ++iter) {
    r=*iter;
    if(disjoint(r,bs)) {
      //std::cerr << "disjoint("<<r<<","<<bs<<")\n";
    } else {
      //std::cerr << "adjoining "<<*iter<<"\n";
      gcls.adjoin(*iter);
    }
  }
  return gcls;
}


template<class R, class BS>
inline
Geometry::GridCellListSet<R>
inner_approximation_of_basic_set(const BS& bs, const Geometry::Grid<R>& g) 
{
  Geometry::GridCellListSet<R> gcls(g);
  Geometry::GridBlock<R> gbb=outer_approximation(bs.bounding_box(),g);
  Geometry::Rectangle<R> r(bs.dimension());
  
  for(typename Geometry::GridBlock<R>::const_iterator iter=gbb.begin(); iter!=gbb.end(); ++iter) {
    r=*iter;
    if(subset(r,bs)) {
      gcls.adjoin(*iter);
    }
  }
  return gcls;
}



}





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
Geometry::GridCellListSet<R>::_instantiate_geometry_operators()
{
  typedef Numeric::Interval<R> I;
  tribool tb;
  //Rectangle<R>* r=0;
  Zonotope<R,R>* z=0;
  Zonotope<I,R>* ez=0;
  Zonotope<I,I>* iz=0;
  Polytope<R>* pltp=0;
  Polyhedron<R>* plhd=0;
  SetInterface<R>* set=0;
  
  Grid<R>* g=0;
  GridBlock<R>* gb=0;
  GridCellListSet<R>* gcls=0;
  
  tb=Geometry::subset(*gcls,*gb);
  
  Geometry::outer_approximation(*pltp,*g);
  Geometry::outer_approximation(*plhd,*g);
  Geometry::outer_approximation(*z,*g);
  Geometry::outer_approximation(*ez,*g);
  Geometry::outer_approximation(*iz,*g);
  Geometry::outer_approximation(*set,*g);
  
  Geometry::inner_approximation(*pltp,*g);
  Geometry::inner_approximation(*plhd,*g);
  Geometry::inner_approximation(*z,*g);
  Geometry::inner_approximation(*ez,*g);
  Geometry::inner_approximation(*iz,*g);
  Geometry::inner_approximation(*set,*g);
  
  Geometry::lower_approximation(*set,*g);
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
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const Polyhedron<R>& ph, const Grid<R>& g) 
{
  return ::outer_approximation_of_basic_set(ph,g);
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const Polyhedron<R>& ph, const Grid<R>& g) 
{
  return ::inner_approximation_of_basic_set(ph,g);
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const Polytope<R>& p, const Grid<R>& g) 
{
  return ::outer_approximation_of_basic_set(p,g);
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const Polytope<R>& p, const Grid<R>& g) 
{
  return ::inner_approximation_of_basic_set(p,g);
}


template<class R,class R0,class R1>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const Zonotope<R0,R1>& z, const Grid<R>& g) 
{
  return ::outer_approximation_of_basic_set(z,g);
}

template<class R,class R0,class R1>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const Zonotope<R0,R1>& z, const Grid<R>& g) 
{
  return ::inner_approximation_of_basic_set(z,g);
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const SetInterface<R>& set, const Grid<R>& g) 
{
  ARIADNE_LOG(4,"GridCellListSet outer_approximation(SetInterface, Grid)\n");
  GridCellListSet<R> result(g);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(set,g,"outer_approximation(SetInterface<R>,Grid<R>)");
  ARIADNE_CHECK_BOUNDED(set,"outer_approximation(SetInterface<R>,Grid<R>)");
  
  const GridBlock<R> gb=outer_approximation(set.bounding_box(),g);
  Rectangle<R> r(g.dimension());
  for(typename GridBlock<R>::const_iterator iter=gb.begin(); iter!=gb.end(); ++iter) {
    r=*iter;
    if(!bool(set.disjoint(r))) {
      result.adjoin(*iter);
    }
  }
  return result;
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const SetInterface<R>& set, const Grid<R>& g) 
{
  ARIADNE_LOG(4,"GridCellListSet inner_approximation(SetInterface, Grid)\n");
  GridCellListSet<R> result(g);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(set,g,"inner_approximation(SetInterface<R>,Grid<R>)\n");
  ARIADNE_CHECK_BOUNDED(set,"inner_approximation(SetInterface<R>,Grid<R>)");
  
  const GridBlock<R> gb=outer_approximation(set.bounding_box(),g);
  Rectangle<R> r(g.dimension());
  for(typename GridBlock<R>::const_iterator iter=gb.begin(); iter!=gb.end(); ++iter) {
    r=*iter;
    if((set.superset(r))) {
      result.adjoin(*iter);
    }
  }
  return result;
}


template<class R>
Geometry::ListSet< Geometry::Rectangle<R> >
Geometry::lower_approximation(const SetInterface<R>& s, const Grid<R>& g) 
{
  ARIADNE_LOG(4,"ListSet<Rectangle> lower_approximation(SetInterface s, Grid fg)\n");
  FiniteGrid<R> fg(g,s.bounding_box());
  return lower_approximation(s,fg);
}











}
