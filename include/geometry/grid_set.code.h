/***************************************************************************
 *            grid_set.code.h
 *
 *  10 January 2005
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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

#include "grid_set.h"

#include <ostream>

#include "../base/stlio.h"

#include "../combinatoric/array_operations.h"

#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/list_set.h"
#include "../geometry/partition_tree_set.h"

#include "../geometry/set_interface.h"

namespace Ariadne {
  namespace Geometry {

    extern int verbosity;
    
    // GridCell ----------------------------------------------------------------
    
    template<class R>
    GridCell<R>::GridCell(const Grid<R>& g, const Combinatoric::LatticeCell& pos)
      : _grid_ptr(&g), _lattice_set(pos)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(g,pos,"GridCell<R>::GridCell(Grid<R>,LatticeCell");
    }

    template<class R>
    GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
      : _grid_ptr(&g), _lattice_set(pos)
    {
      ARIADNE_CHECK_DIMENSION(g,pos.size(),"GridCell::GridCell(Grid,IndexArray)");
    }

    
    template<class R>
    R
    GridCell<R>::lower_bound(dimension_type i) const 
    {
      return _grid_ptr->subdivision_coordinate(i,_lattice_set.lower_bound(i));
    }
    
    
    template<class R>
    R
    GridCell<R>::upper_bound(dimension_type i) const 
    {
      return _grid_ptr->subdivision_coordinate(i,_lattice_set.upper_bound(i));
    }
    
    template<class R>
    GridBlock<R>
    GridCell<R>::neighbourhood() const 
    {
      return GridBlock<R>(this->_grid_ptr,this->_lattice_set.neighbourhood());
    }
    

    template<class R>
    void
    GridCell<R>::_instantiate_geometry_operators() 
    {    
    }
    
    
    
    
    // GridBlock --------------------------------------------------------------
    
    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    { 
      _lattice_set.set_lower_bound(0,1);
      //_lattice_set.set_lower_bound(0,0);
      _lattice_set.set_upper_bound(0,0);
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const Combinatoric::LatticeBlock& lb)
      : _grid_ptr(&g), _lattice_set(lb)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(g,lb,"GridBlock::GridBlock(Grid g, LatticeBlock lb)");
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid_ptr(&g), _lattice_set(l,u)
    {
      ARIADNE_CHECK_DIMENSION(g,l.size(),"GridBlock::GridBlock(Grid g, IndexArray l, IndexArray u)");
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const Rectangle<R>& r)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(g,r,"GridBlock::GridBlock(Grid g,Rectangle r)");
      for(dimension_type i=0; i!=dimension(); ++i) {
        /* TODO: Catch and rethrow exceptions */
        _lattice_set.set_lower_bound(i,g.subdivision_index(i,r.lower_bound(i)));
        _lattice_set.set_upper_bound(i,g.subdivision_index(i,r.upper_bound(i)));
      }
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const GridCell<R>& gc)
      : _grid_ptr(gc._grid_ptr), _lattice_set(gc.lattice_set())
    {
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const GridBlock<R>& gb)
      : _grid_ptr(gb._grid_ptr), _lattice_set(gb.lattice_set())
    {
    }
    
    
    template<class R>
    GridBlock<R>&
    GridBlock<R>::operator=(const GridBlock<R>& gb)
    {
      if(this!=&gb) {
        this->_grid_ptr = gb._grid_ptr;
        this->_lattice_set=gb._lattice_set;
      }
      return *this;
    }
    
    
    template<class R>
    R
    GridBlock<R>::lower_bound(dimension_type i) const 
    {
      return _grid_ptr->subdivision_coordinate(i,_lattice_set.lower_bound(i));
    }
    
    
    template<class R>
    R
    GridBlock<R>::upper_bound(dimension_type i) const 
    {
      return _grid_ptr->subdivision_coordinate(i,_lattice_set.upper_bound(i));
    }
    
    template<class R>
    GridBlock<R>
    GridBlock<R>::neighbourhood() const 
    {
      return GridBlock<R>(this->_grid_ptr,this->_lattice_set.neighbourhood());
    }
    
    
    template<class R>
    void
    GridBlock<R>::_instantiate_geometry_operators() 
    {
      tribool tb;
      Rectangle<R>* r=0;
      Grid<R>* g=0;
      GridCell<R>* gc=0;
      GridBlock<R>* gb=0;
      
      tb=Geometry::subset(*r,*gb);
      
      tb=Geometry::overlap(*gb,*gb);
      tb=Geometry::subset(*gc,*gb);
      tb=Geometry::subset(*gb,*gb);

      *gb=Geometry::under_approximation(*r,*g);
    }
    
    
    
    
    // GridCellListSet ------------------------------------------------------------

    template<class R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
      : _grid_ptr(new Grid<R>(g)), _lattice_set(g.dimension())
    {
    }

    
    template<class R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g, 
                                        const Combinatoric::LatticeCellListSet& lcls)
      : _grid_ptr(new Grid<R>(g)), _lattice_set(lcls)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(g,lcls,"GridCellListSet::GridCellListSet(Grid g, LatticeCellListSet lcls)");
    }

    
    template<class R>
    GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
      : _grid_ptr(gms._grid_ptr), _lattice_set(gms.dimension())
    {
      this->_lattice_set.adjoin(gms._lattice_set);
    }

    
    template<class R>
    GridCellListSet<R>::GridCellListSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(gcls._grid_ptr), _lattice_set(gcls._lattice_set)
    {
    }

    template<class R>
    GridCellListSet<R>&
    GridCellListSet<R>::operator=(const GridCellListSet<R>& gcls)
    {
      if(this!=&gcls) {
        this->_grid_ptr = gcls._grid_ptr;
        this->_lattice_set=gcls._lattice_set;
      }
      return *this;
    }

    

    

    template<class R>
    GridCellListSet<R>::operator ListSet< Rectangle<R> >() const
    {
      ListSet< Rectangle<R> > result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }

    
    template<class R>
    void
    GridCellListSet<R>::clear()
    {
      this->_lattice_set.clear();
    }
    
    
    template<class R>
    void
    GridCellListSet<R>::_instantiate_geometry_operators()
    {
      tribool tb;
      Zonotope<R>* z=0;
      Zonotope< Numeric::Interval<R> >* iz=0;
      Polytope<R>* ply=0;
      Polyhedron<R>* phd=0;
      Grid<R>* g=0;
      GridBlock<R>* gb=0;
      GridCellListSet<R>* gcls=0;

      tb=Geometry::subset(*gcls,*gb);
      
      *gcls=Geometry::over_approximation(*z,*g);
      *gcls=Geometry::over_approximation(*ply,*g);
      *gcls=Geometry::over_approximation(*phd,*g);
      *gcls=Geometry::over_approximation(*iz,*g);

      *gcls=Geometry::under_approximation(*z,*g);
      *gcls=Geometry::under_approximation(*ply,*g);
      *gcls=Geometry::under_approximation(*phd,*g);
    }


    // GridMaskSet ------------------------------------------------------------
    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg)
      : _grid_ptr(new Grid<R>(fg.grid())), _lattice_set(fg.lattice_block()) 
    { 
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg, const BooleanArray& m)
      : _grid_ptr(new Grid<R>(fg.grid())), _lattice_set(fg.lattice_block(),m)
    {
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Rectangle<R>& bb)
      : _grid_ptr(new Grid<R>(g)), _lattice_set(g.index_block(bb)) 
    { 
    }

    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b)
      : _grid_ptr(new Grid<R>(g)), _lattice_set(b) 
    { 
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b, const BooleanArray& m)
      : _grid_ptr(new Grid<R>(g)), _lattice_set(b,m)
    {
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeMaskSet& ms)
      : _grid_ptr(new Grid<R>(g)), _lattice_set(ms)
    {
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
      : _grid_ptr(gms._grid_ptr), _lattice_set(gms._lattice_set)
    {
    }
    
    
    template<class R>
    GridMaskSet<R>&
    GridMaskSet<R>::operator=(const GridMaskSet<R>& gms) 
    {
      if(this!=&gms) {
        this->_grid_ptr=gms._grid_ptr;
        this->_lattice_set=gms._lattice_set;
      }
      return *this;
    }
    
    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(gcls._grid_ptr),
        _lattice_set(gcls.lattice_set())
    {
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
      return !Geometry::disjoint(*this,Rectangle<R>(pt));
    }


    template<class R>
    tribool
    GridMaskSet<R>::superset(const Rectangle<R>& r) const
    {
      return Geometry::subset(r,*this);
    }


    template<class R>
    tribool
    GridMaskSet<R>::intersects(const Rectangle<R>& r) const
    {
      return !Geometry::disjoint(*this,r);
    }


    template<class R>
    tribool
    GridMaskSet<R>::disjoint(const Rectangle<R>& r) const
    {
      return Geometry::disjoint(*this,r);
    }


    template<class R>
    tribool
    GridMaskSet<R>::subset(const Rectangle<R>& r) const
    {
      return Geometry::subset(*this,r);
    }

    
    template<class R> 
    Rectangle<R> 
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
      tribool tb;
      Point< Numeric::Interval<R> >* ipt=0;
      Rectangle<R>* r=0;
      ListSet< Rectangle<R> >* rls=0;
      ListSet< Zonotope<R> >* zls=0;
      Grid<R>* g=0;
      FiniteGrid<R>* fg=0;
      GridCell<R>* gc=0;
      GridBlock<R>* gb=0;
      GridCellListSet<R>* gcls=0;
      GridMaskSet<R>* gms=0;
      PartitionTreeSet<R>* pts=0;
      SetInterface<R>* set=0;
      
      tb=Geometry::subset(*r,*gms);
      tb=Geometry::subset(*gms,*r);
      tb=Geometry::disjoint(*r,*gms);
      tb=Geometry::disjoint(*gms,*r);
       
      tb=Geometry::overlap(*gb,*gms);
      tb=Geometry::overlap(*gms,*gb);
      tb=Geometry::overlap(*gms,*gms);
    
      tb=Geometry::subset(*gc,*gms);
      tb=Geometry::subset(*gb,*gms);
      tb=Geometry::subset(*gcls,*gms);
      tb=Geometry::subset(*gms,*gms);
      
      *gms=Geometry::regular_intersection(*gb,*gms);
      *gms=Geometry::regular_intersection(*gms,*gb);
      *gcls=Geometry::regular_intersection(*gcls,*gms);
      *gcls=Geometry::regular_intersection(*gms,*gcls);
      *gms=Geometry::regular_intersection(*gms,*gms);
      *gcls=Geometry::difference(*gcls,*gms);
      *gms=Geometry::difference(*gms,*gms);
      *gms=Geometry::join(*gms,*gms);

      *gb=Geometry::over_approximation(*ipt,*g);
  
      *gms=Geometry::over_approximation(*rls,*fg);
      *gms=Geometry::over_approximation(*zls,*fg);
      *gms=Geometry::over_approximation(*gms,*fg);
      *gms=Geometry::over_approximation(*pts,*fg);
      *gms=Geometry::over_approximation(*set,*fg);
      
      *gms=Geometry::over_approximation(*set,*g);

      *gms=Geometry::under_approximation(*rls,*fg);
      *gms=Geometry::under_approximation(*gms,*fg);
      *gms=Geometry::under_approximation(*pts,*fg);
      *gms=Geometry::under_approximation(*set,*fg);
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
    disjoint(const Rectangle<R>& r, const GridMaskSet<R>& gms) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gms,"tribool disjoint(Rectangle r, GridMaskSet gms)");
      Rectangle<R> br=closed_intersection(r,Rectangle<R>(gms.bounding_box()));
      GridBlock<R> gb=outer_approximation(br,gms.grid());
      return !overlap(gb,gms);
    }
    

    template<class R>
    tribool
    disjoint(const GridMaskSet<R>& gms, const Rectangle<R>& r) {
      return disjoint(r,gms);
    }
    
    

    template<class R>
    tribool
    overlap(const GridBlock<R>& gb1, const GridBlock<R>& gb2) {
      ARIADNE_CHECK_SAME_GRID(gb1,gb2,"tribool overlap(GridBlock gb1,GridBlock gb2)");
      return overlap(gb1.lattice_set(),gb2.lattice_set());
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
    overlap(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
    {
      ARIADNE_CHECK_SAME_GRID(gms1,gms2,"tribool overlap(GridMaskSet gms2, GridMaskSet gms2)");
      return overlap(gms1.lattice_set(),gms2.lattice_set());
    }




    template<class R>
    tribool
    subset(const GridCell<R>& gc, const GridBlock<R>& gb)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(gc,gb,"tribool subset(GridCell gc, GridBlock gb)");
      if(gc.grid()==gb.grid()) {
        return subset(gc.lattice_set(),gb.lattice_set());
      }
      return subset(Rectangle<R>(gc),Rectangle<R>(gb));
    }

    template<class R>
    tribool
    subset(const GridBlock<R>& gb1, const GridBlock<R>& gb2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(gb1,gb2,"tribool subset(GridBlock gb1, GridBlock gb2)");
      if(gb1.grid()==gb2.grid()) {
        return subset(gb1.lattice_set(),gb2.lattice_set());
      }
      return subset(Rectangle<R>(gb1),Rectangle<R>(gb2));
    }

    template<class R>
    tribool
    subset(const GridCellListSet<R>& gcls, const GridBlock<R>& gb)
    {
      ARIADNE_CHECK_SAME_GRID(gcls,gb,"tribool subset(GridCellListSet gcls, GridBlock gb)");
      return subset(gcls.lattice_set(),gb.lattice_set());
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
    subset(const Rectangle<R>& r, const GridBlock<R>& gb)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gb,"tribool subset(Rectangle r, GridBlock gb)");
      return subset(r,Rectangle<R>(gb));
    }
    
    

    template<class R>
    tribool
    subset(const Rectangle<R>& r, const GridMaskSet<R>& gms)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gms,"tribool subset(Rectangle r, GridMaskSet gms)");
      if(!subset(r,gms.bounding_box())) {
        return false;
      }
      return subset(over_approximation(r,gms.grid()),gms);
    }


    template<class R>
    tribool
    subset(const GridMaskSet<R>& gms, const Rectangle<R>& r)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(gms,r,"tribool subset(GridMaskSet gms, Rectangle r)");
      if(!subset(gms,r.bounding_box())) {
        return false;
      }
      return subset(gms,over_approximation(r,gms.grid()));
    }





    template<class R>
    GridMaskSet<R>
    join(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
    {
      ARIADNE_CHECK_SAME_GRID(gms1,gms2,"GridMaskSet join(GridMaskSet gms1, GridMaskSet gms2)");
      if(gms1.block()==gms2.block()) { 
        return GridMaskSet<R>(gms1.grid(), join(gms1.lattice_set(),gms2.lattice_set()));
      } else {
        GridMaskSet<R> result(gms1.grid(),Combinatoric::rectangular_hull(gms1.block(),gms2.block()));
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
        GridMaskSet<R> result(gms1.grid(),Combinatoric::rectangular_hull(gms1.block(),gms2.block()));
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
        GridMaskSet<R> result(gms1.grid(),Combinatoric::rectangular_hull(gms1.block(),gms2.block()));
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
      return GridCellListSet<R>(gcls.grid(),Combinatoric::difference(gcls.lattice_set(),gms.lattice_set()));
    }

    template<class R>
    GridMaskSet<R>::operator ListSet< Rectangle<R> >() const
    {
      ARIADNE_CHECK_BOUNDED(*this,"GridMaskSet::operator ListSet<Rectangle>()");
      ListSet< Rectangle<R> > result(this->dimension());
      for(typename GridMaskSet::const_iterator riter=begin(); riter!=end(); ++riter) {
        Rectangle<R> r(*riter);
        result.push_back(r);
      }
      return result;
    }





  
  

    // Over and under approximations ------------------------------------------
    
    template<class R>
    GridBlock<R>
    outer_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      if(r.empty()) {
        return GridBlock<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_upper_index(i,r.lower_bound(i))-1;
        upper[i]=g.subdivision_lower_index(i,r.upper_bound(i))+1;
      }
      return GridBlock<R>(g,lower,upper);
    }
    
    template<class R>
    GridBlock<R>
    inner_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      if(r.empty()) {
        return GridBlock<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,r.lower_bound(i))+1;
        upper[i]=g.subdivision_upper_index(i,r.upper_bound(i))-1;
      }
      return GridBlock<R>(g,lower,upper);
    }
    
    template<class R>
    GridBlock<R>
    over_approximation(const Point< Numeric::Interval<R> >& ipt, const Grid<R>& g) 
    {
      return over_approximation(Rectangle<R>(ipt),g);
    }
    
    
    
    template<class R>
    GridBlock<R>
    over_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      if(r.empty()) {
        return GridBlock<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,r.lower_bound(i));
        upper[i]=g.subdivision_upper_index(i,r.upper_bound(i));
      }

      return GridBlock<R>(g,lower,upper);
    }
    
    template<class R>
    GridBlock<R>
    under_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      
      if(r.empty()) {
        return GridBlock<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,r.lower_bound(i))+1;
        upper[i]=g.subdivision_upper_index(i,r.upper_bound(i))-1;
      }

      return GridBlock<R>(g,lower,upper);
    }
    
    template<class R>
    GridCellListSet<R>
    over_approximation(const Zonotope<R>& z, const Grid<R>& g) 
    {
      if(verbosity>7) { std::clog << "GridCellListSet over_approximation(Zonotope z, Grid g)" << std::endl; }
      GridCellListSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(z,g,"over_approximation(Zonotope<R>,Grid<R>)");
      if(z.empty()) {
        return result; 
      }
      Rectangle<R> bb=z.bounding_box();

      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(disjoint(Rectangle<R>(cell),z)) {
        } else {
          result.adjoin(cell);
        }
      }

      return result;
    }


    template<class R>
    GridCellListSet<R>
    over_approximation(const Polytope<R>& p, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(p,g,"over_approximation(Polytope<R>,Grid<R>)");
      if(p.empty()) {
        return result; 
      }
      Rectangle<R> bb=p.bounding_box();

      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(disjoint(Rectangle<R>(cell),p)) {
        } else {
          result.adjoin(cell);
        }
      }

      return result;
    }

   
    template<class R>
    GridCellListSet<R>
    over_approximation(const Polyhedron<R>& p, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(p,g,"over_approximation(Polyhedron<R>,Grid<R>)");
      
      // omit emptyness and boundedness checks
          
      Rectangle<R> bb=p.bounding_box();

      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(disjoint(Rectangle<R>(cell),p)) {
        } else {
          result.adjoin(cell);
        }
      }

      return result;
    }


    template<class R>
    GridCellListSet<R>
    over_approximation(const Zonotope< Numeric::Interval<R> >& z, const Grid<R>& g) 
    {
      if(verbosity>7) { std::clog << "GridCellListSet over_approximation(Zonotope<Interval>, Grid g)" << std::endl; }
      GridCellListSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(z,g,"over_approximation(Zonotope<R>,Grid<R>)");
      if(z.empty()) {
        return result; 
      }
      
      size_type ng=z.number_of_generators();
      Point< Numeric::Interval<R> > c=z.centre();
      LinearAlgebra::Matrix< Numeric::Interval<R> > G=z.generators();
      
      Point<R> ac=approximate_value(c);
      LinearAlgebra::Matrix<R> aG=approximate_value(G);
      LinearAlgebra::Vector< Numeric::Interval<R> > e=(c-ac)+(G-aG)*LinearAlgebra::Vector< Numeric::Interval<R> >(ng,Numeric::Interval<R>(-1,1));
      
      Zonotope<R> az(ac,aG);
      
     
      Rectangle<R> bb=z.bounding_box();

      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(disjoint(Rectangle<R>(cell)+e,az)) {
        } else {
          result.adjoin(cell);
        }
      }

      return result;
    }

    template<class R>
    GridCellListSet<R>
    under_approximation(const Zonotope<R>& z, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(z,g,"under_approximation(Zonotope<R>,Grid<R>)");
      if(z.empty()) {
        return result; 
      }
      
      Rectangle<R> bb=z.bounding_box();
      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(subset(Rectangle<R>(cell),z)) {
          result.adjoin(cell);
        }
      }
      return result;
    }

    template<class R>
    GridCellListSet<R>
    under_approximation(const Polytope<R>& p, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(p,g,"under_approximation(Polytope<R>,Grid<R>)");
      if(p.empty()) {
        return result; 
      }
      
      Rectangle<R> bb=p.bounding_box();
      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(subset(Rectangle<R>(cell),p)) {
          result.adjoin(cell);
        }
      }
      return result;
    }

    template<class R>
    GridCellListSet<R>
    under_approximation(const Polyhedron<R>& p, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(p,g,"under_approximation(Polyhedron<R>,Grid<R>)");
      
      Rectangle<R> bb=p.bounding_box();
      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(subset(Rectangle<R>(cell),p)) {
          result.adjoin(cell);
        }
      }
      return result;
    }

    
    
    template<class R>
    GridMaskSet<R>
    over_approximation(const GridMaskSet<R>& gms, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(gms,g,"over_approximation(GridMaskSet<R>,FiniteGrid<R>)");

      for(typename GridMaskSet<R>::const_iterator iter=gms.begin(); iter!=gms.end(); ++iter) {
        result.adjoin(over_approximation(Rectangle<R>(*iter),g.grid()));
      }
      return result;
    }
    
    template<class R>
    GridMaskSet<R>
    over_approximation(const PartitionTreeSet<R>& pts, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(pts,g,"over_approximation(PartitionTreeSet<R>,FiniteGrid<R>)");

      for(typename PartitionTreeSet<R>::const_iterator iter=pts.begin(); iter!=pts.end(); ++iter) {
        result.adjoin(over_approximation(Rectangle<R>(*iter),g.grid()));
      }
      return result;
    }
    
    template<class R>
    GridMaskSet<R>
    over_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg) 
    {
      GridMaskSet<R> result(fg);
      ARIADNE_CHECK_EQUAL_DIMENSIONS(set,fg,"over_approximation(PartitionTreeSet<R>,FiniteGrid<R>)");
      
      const Grid<R>& g=fg.grid();
      const Combinatoric::LatticeBlock& lb=fg.lattice_block();
      
      for(typename Combinatoric::LatticeBlock::const_iterator iter=lb.begin(); iter!=lb.end(); ++iter) {
        GridCell<R> gc(g,*iter);
        if(!bool(set.disjoint(Rectangle<R>(gc)))) {
          result.adjoin(gc);
        }
      }
      return result;
    }
    
    
    template<class R>
    GridMaskSet<R>
    over_approximation(const SetInterface<R>& set, const Grid<R>& g) 
    {
      FiniteGrid<R> fg(g,set.bounding_box());
      return over_approximation(set,fg);
    }
      

    template<class R>
    GridMaskSet<R>
    under_approximation(const GridMaskSet<R>& gms, const FiniteGrid<R>& g) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(gms,g,"under_approximation(GridMaskSet<R>,FiniteGrid<R>)");
      
      GridMaskSet<R> result(g);
      GridMaskSet<R> bb(g);
      bb.adjoin(bb.bounds());
      Rectangle<R> r;
      
      for(typename GridMaskSet<R>::const_iterator iter=bb.begin(); iter!=bb.end(); ++iter) {
        r=*iter;
        if(subset(r,gms)) {
          const GridCell<R>& cell=*iter;
          result.adjoin(cell);
        }
      }
      return result;
    }

    template<class R>
    GridMaskSet<R>
    under_approximation(const PartitionTreeSet<R>& pts, const FiniteGrid<R>& g) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    GridMaskSet<R>
    under_approximation(const SetInterface<R>& s, const FiniteGrid<R>& fg) 
    {
      GridMaskSet<R> result(fg);
      GridBlock<R> gb(fg.grid(),fg.lattice_block());
      Rectangle<R> r;
      for(typename GridBlock<R>::const_iterator gb_iter=gb.begin();
          gb_iter!=gb.end(); ++gb_iter)
      {
        r=*gb_iter;
        if(s.superset(r)) {
          result.adjoin(fg);
        }
      }
      return result;
    }
    

    
    template<class R>
    inline
    GridBlock<R>
    over_approximation(const Rectangle<R>& r, const GridMaskSet<R>& gms) 
    {
      return over_approximation(r,gms.grid());
    }
    
    template<class R>
    inline
    GridBlock<R>
    under_approximation(const Rectangle<R>& r, const GridMaskSet<R>& gms) 
    {
      return under_approximation(r,gms.grid());
    }
    



    // Input/output------------------------------------------------------------

    template<class R>
    std::ostream&
    GridCell<R>::write(std::ostream& os) const 
    {
      return os << Rectangle<R>(*this);
    }

    
    template<class R>
    std::ostream&
    GridBlock<R>::write(std::ostream& os) const 
    {
      return os << Rectangle<R>(*this);
    }

    
    template<class R>
    std::ostream&     
    GridCellListSet<R>::write(std::ostream& os) const 
    {
      os << "GridCellListSet<" << Numeric::name<R>() << ">(\n";
      os << "  grid=" << this->grid() << ",\n";
      os << "  size=" << this->size() << ",\n";
      os << "  lattice_set=" << this->lattice_set();
      os << ")" << std::endl;
      return os;
    }



    template<class R>
    std::ostream& 
    GridMaskSet<R>::write(std::ostream& os) const 
    {
      os << "GridMaskSet<" << Numeric::name<R>() << ">("<< std::endl;
      os << "  grid=" << this->grid() << "," << std::endl;
      os << "  bounds=" << this->bounds() << "," << std::endl;
      os << "  size=" << this->size() << "," << std::endl;
      os << "  capacity=" << this->capacity() << "," << std::endl;
      os << "  mask=" << this->mask() << std::endl;
      os << ")\n";
      return os;
    }

  }
}
