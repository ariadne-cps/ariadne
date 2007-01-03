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
#include "../geometry/parallelotope.h"
#include "../geometry/polytope.h"
#include "../geometry/list_set.h"
#include "../geometry/partition_tree_set.h"
#include "../geometry/set.h"

namespace Ariadne {
  namespace Geometry {

    extern int verbosity;
    
    // GridCell ----------------------------------------------------------------
    
    template<class R>
    GridCell<R>::GridCell(const Grid<R>& g, const Combinatoric::LatticeCell& pos)
      : _grid(g), _lattice_set(pos)
    {
      check_equal_dimensions(g,pos,"GridCell<R>::GridCell(Grid<R>,LatticeCell");
    }

    template<class R>
    GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
      : _grid(g), _lattice_set(pos)
    {
      check_dimension(g,pos.size(),__PRETTY_FUNCTION__);
    }

    
    template<class R>
    R
    GridCell<R>::lower_bound(dimension_type i) const 
    {
      return _grid.subdivision_coordinate(i,_lattice_set.lower_bound(i));
    }
    
    
    template<class R>
    R
    GridCell<R>::upper_bound(dimension_type i) const 
    {
      return _grid.subdivision_coordinate(i,_lattice_set.upper_bound(i));
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
    GridBlock<R>::GridBlock(const Grid<R>& g, const Combinatoric::LatticeBlock& b)
      : _grid_ptr(&g), _lattice_set(b)
    {
      check_equal_dimensions(g,b,"GridBlock<R>::GridBlock(Grid<R>,LatticeBlock)");
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid_ptr(&g), _lattice_set(l,u)
    {
      check_dimension(g,l.size(),__PRETTY_FUNCTION__);
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const Rectangle<R>& r)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
      check_equal_dimensions(g,r,"GridBlock<R>::GridBlock(Grid<R>,Rectangle<R>)");
      for(dimension_type i=0; i!=dimension(); ++i) {
        /* TODO: Catch and rethrow exceptions */
        _lattice_set.set_lower_bound(i,g.subdivision_index(i,r.lower_bound(i)));
        _lattice_set.set_upper_bound(i,g.subdivision_index(i,r.upper_bound(i)));
      }
    }
    
    
    template<class R>
    GridBlock<R>::GridBlock(const GridCell<R>& gc)
      : _grid_ptr(&gc.grid()), _lattice_set(gc.lattice_set())
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
        this->_grid_ptr=gb._grid_ptr;
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
    void
    GridBlock<R>::_instantiate_geometry_operators() 
    {
      tribool tb;
      Rectangle<R>* r=0;
      Grid<R>* g=0;
      GridCell<R>* gc=0;
      GridBlock<R>* gb=0;
      
      tb=subset(*r,*gb);
      
      tb=overlap(*gb,*gb);
      tb=subset(*gc,*gb);
      tb=subset(*gb,*gb);

      *gb=under_approximation(*r,*g);
    }
    
    
    
    
    // GridCellListSet ------------------------------------------------------------

    template<class R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
    }

    
    template<class R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g, 
                                        const Combinatoric::LatticeCellListSet& lcls)
      : _grid_ptr(&g), _lattice_set(lcls)
    {
      check_equal_dimensions(g,lcls,"GridCellListSet<R>::GridCellListSet(Grid<R>,LatticeCellListSet)");
    }

    
    template<class R>
    GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
      : _grid_ptr(&gms.grid()), _lattice_set(gms.dimension())
    {
      this->_lattice_set.adjoin(gms._lattice_set);
    }

    
    template<class R>
    GridCellListSet<R>::GridCellListSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()), _lattice_set(gcls._lattice_set)
    {
    }

    
    // FIXME: Memory leak
    template<class R>
    GridCellListSet<R>::GridCellListSet(const ListSet<R,Rectangle>& rls)
      : _grid_ptr(new IrregularGrid<R>(rls)), _lattice_set(rls.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator list_set_const_iterator;
      for(list_set_const_iterator iter=rls.begin(); iter!=rls.end(); ++iter) {
        this->adjoin(GridBlock<R>(grid(),*iter));
      }
    }

    
    template<class R>
    GridCellListSet<R>::operator ListSet<R,Rectangle>() const
    {
      ListSet<R,Rectangle> result(dimension());
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
      Parallelotope<R>* pl=0;
      Polytope<R>* p=0;
      Polyhedron<R>* ph=0;
      Grid<R>* g=0;
      GridBlock<R>* gb=0;
      GridCellListSet<R>* gcls=0;

      tb=subset(*gcls,*gb);
      
      *gcls=over_approximation(*z,*g);
      *gcls=over_approximation(*pl,*g);
      *gcls=over_approximation(*p,*g);

      *gcls=under_approximation(*z,*g);
      *gcls=under_approximation(*pl,*g);
      *gcls=under_approximation(*p,*g);
      *gcls=under_approximation(*ph,*g);
    }


    // GridMaskSet ------------------------------------------------------------
    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& g)
      : _grid_ptr(&g.grid()), _lattice_set(g.bounds()) 
    { 
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg, const BooleanArray& m)
      : _grid_ptr(&fg.grid()), _lattice_set(fg.bounds(),m)
    {
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b)
      : _grid_ptr(&g), _lattice_set(b) 
    { 
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b, const BooleanArray& m)
      : _grid_ptr(&g), _lattice_set(b,m)
    {
      const IrregularGrid<R>* irregular_grid_ptr=dynamic_cast<const IrregularGrid<R>*>(_grid_ptr);
      if(irregular_grid_ptr) {
        if(!subset(b,irregular_grid_ptr->block())) {
          throw std::runtime_error("Lattice block does not lie in grid bounds");
        }
      }
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeMaskSet& ms)
      : _grid_ptr(&g), _lattice_set(ms)
    {
      const IrregularGrid<R>* irregular_grid_ptr=dynamic_cast<const IrregularGrid<R>*>(_grid_ptr);
      if(irregular_grid_ptr) {
        if(!subset(ms.block(),irregular_grid_ptr->block())) {
          throw std::runtime_error("Lattice block does not lie in grid bounds");
        }
      }
    }

    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
      : _grid_ptr(&gms.grid()), _lattice_set(gms._lattice_set)
    {
    }
    
    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()),
        _lattice_set(gcls.lattice_set())
    {
    }
    
    
    // FIXME: Memory leak
    template<class R>
    GridMaskSet<R>::GridMaskSet(const ListSet<R,Rectangle>& rls) 
      : _grid_ptr(new IrregularGrid<R>(rls)), 
        _lattice_set(dynamic_cast<const IrregularGrid<R>*>(_grid_ptr)->block())
    {
      //FIXME: Memory leak!    

      for(typename ListSet<R,Rectangle>::const_iterator riter=rls.begin(); 
          riter!=rls.end(); ++riter) 
      {
        GridBlock<R> r(grid(),*riter);
        adjoin(r);
      }
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
      Rectangle<R>* r=0;
      ListSet<R,Rectangle>* rls=0;
      ListSet<R,Zonotope>* zls=0;
      ListSet<R,Parallelotope>* plls=0;
      Grid<R>* g=0;
      FiniteGrid<R>* fg=0;
      GridCell<R>* gc=0;
      GridBlock<R>* gb=0;
      GridCellListSet<R>* gcls=0;
      GridMaskSet<R>* gms=0;
      PartitionTreeSet<R>* pts=0;
      Set<R>* set=0;
      
      tb=subset(*r,*gms);
      tb=disjoint(*r,*gms);
      tb=disjoint(*gms,*r);
       
      tb=overlap(*gb,*gms);
      tb=overlap(*gms,*gb);
      tb=overlap(*gms,*gms);
    
      tb=subset(*gc,*gms);
      tb=subset(*gb,*gms);
      tb=subset(*gcls,*gms);
      tb=subset(*gms,*gms);

      *gms=regular_intersection(*gb,*gms);
      *gms=regular_intersection(*gms,*gb);
      *gcls=regular_intersection(*gcls,*gms);
      *gcls=regular_intersection(*gms,*gcls);
      *gms=regular_intersection(*gms,*gms);
      *gcls=difference(*gcls,*gms);
      *gms=difference(*gms,*gms);
      *gms=join(*gms,*gms);

      *gms=over_approximation(*rls,*fg);
      *gms=over_approximation(*zls,*fg);
      *gms=over_approximation(*plls,*fg);
      *gms=over_approximation(*gms,*fg);
      *gms=over_approximation(*pts,*fg);
      *gms=over_approximation(*set,*fg);
      
      *gms=over_approximation(*set,*g);

      *gms=under_approximation(*rls,*fg);
      *gms=under_approximation(*gms,*fg);
      *gms=under_approximation(*pts,*fg);
    }
    
    
    
    
    // Geometric predicates ---------------------------------------------------

    template<class R>
    tribool
    disjoint(const GridBlock<R>& A, const GridBlock<R>& B) {
      check_same_grid(A,B,"disjoint(GridBlock<R>,GridBlock<R>)");
      return disjoint(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    disjoint(const GridBlock<R>& A, const GridMaskSet<R>& B) {
      check_same_grid(A,B,"disjoint(GridBlock<R>,GridMaskSet<R>)");
      return disjoint(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    disjoint(const GridMaskSet<R>& A, const GridBlock<R>& B) {
      check_same_grid(A,B,"disjoint(GridMaskSet<R>,GridBlock<R>)");
      return disjoint(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    disjoint(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"disjoint(GridMaskSet<R>,GridMaskSet<R>)");
      return disjoint(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    disjoint(const Rectangle<R>& A, const GridMaskSet<R>& B) {
      check_equal_dimensions(A,B,"disjoint(Rectangle<R>,GridMaskSet<R>)");
      Rectangle<R> r=closed_intersection(A,Rectangle<R>(B.bounding_box()));
      GridBlock<R> gr=outer_approximation(r,B.grid());
      return !overlap(gr,B);
    }
    

    template<class R>
    tribool
    disjoint(const GridMaskSet<R>& A, const Rectangle<R>& B) {
      return disjoint(B,A);
    }
    
    

    template<class R>
    tribool
    overlap(const GridBlock<R>& A, const GridBlock<R>& B) {
      check_same_grid(A,B,"overlap(GridBlock<R>,GridBlock<R>)");
      return overlap(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    overlap(const GridBlock<R>& A, const GridMaskSet<R>& B) {
      check_same_grid(A,B,"overlap(GridBlock<R>,GridMaskSet<R>)");
      return overlap(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    overlap(const GridMaskSet<R>& A, const GridBlock<R>& B) {
      check_same_grid(A,B,"overlap(GridMaskSet<R>,GridBlock<R>)");
      return overlap(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    overlap(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"overlap(GridMaskSet<R>,GridMaskSet<R>)");
      return overlap(A.lattice_set(),B.lattice_set());
    }




    template<class R>
    tribool
    subset(const GridCell<R>& A, const GridBlock<R>& B)
    {
      check_equal_dimensions(A,B,"subset(GridCell<R>,GridBlock<R>)");
      if(A.grid()==B.grid()) {
        return subset(A.lattice_set(),B.lattice_set());
      }
      return subset(Rectangle<R>(A),Rectangle<R>(B));
    }

    template<class R>
    tribool
    subset(const GridBlock<R>& A, const GridBlock<R>& B)
    {
      check_equal_dimensions(A,B,"subset(GridBlock<R>,GridBlock<R>)");
      if(A.grid()==B.grid()) {
        return subset(A.lattice_set(),B.lattice_set());
      }
      return subset(Rectangle<R>(A),Rectangle<R>(B));
    }

    template<class R>
    tribool
    subset(const GridCellListSet<R>& A, const GridBlock<R>& B)
    {
      check_same_grid(A,B,"subset(GridCellListSet<R>,GridBlock<R>)");
      return subset(A.lattice_set(),B.lattice_set());
    }

    template<class R>
    tribool
    subset(const GridMaskSet<R>& A, const GridBlock<R>& B)
    {
      check_same_grid(A,B,"subset(GridMaskSet<R>,GridBlock<R>)");
      return subset(A.lattice_set(),B.lattice_set());
    }

    template<class R>
    tribool
    subset(const GridCell<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"subset(GridCell<R>,GridMaskSet<R>)");
      return subset(A.lattice_set(),B.lattice_set()); 
    }

    template<class R>
    tribool
    subset(const GridBlock<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"subset(GridBlock<R>,GridMaskSet<R>)");
      return subset(A.lattice_set(),B.lattice_set()); 
    }

    template<class R>
    tribool
    subset(const GridCellListSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"subset(GridCellListSet<R>,GridMaskSet<R>)");
      return subset(A.lattice_set(),B.lattice_set()); 
    }

    template<class R>
    tribool
    subset(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"subset(GridMaskSet<R>,GridMaskSet<R>)");
      return subset(A.lattice_set(),B.lattice_set());
    }


    template<class R>
    tribool
    subset(const Rectangle<R>& A, const GridBlock<R>& B)
    {
      check_equal_dimensions(A,B,"subset(Rectangle<R>,GridBlock<R>)");
      return subset(A,Rectangle<R>(B));
    }
    
    

    template<class R>
    tribool
    subset(const Rectangle<R>& A, const GridMaskSet<R>& B)
    {
      check_equal_dimensions(A,B,"subset(Rectangle<R>,GridMaskSet<R>)");
      if(!subset(A,B.bounding_box())) {
        return false;
      }
      return subset(over_approximation(A,B.grid()),B);
    }





    template<class R>
    GridMaskSet<R>
    join(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"join(GridMaskSet<R>,GridMaskSet<R>)");
      if(A.block()!=B.block()) { throw IncompatibleGrids(__PRETTY_FUNCTION__); }
      return GridMaskSet<R>(A.grid(), join(A.lattice_set(),B.lattice_set()));
    }


    template<class R>
    GridMaskSet<R>
    regular_intersection(const GridMaskSet<R>& A, const GridBlock<R>& B)
    {
      check_same_grid(A,B,"regular_intersection(GridMaskSet<R>,GridBlock<R>)");
      return GridMaskSet<R>(A.grid(), regular_intersection(A.lattice_set(),B.lattice_set()));
    }

    template<class R>
    GridMaskSet<R>
    regular_intersection(const GridBlock<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"regular_intersection(GridBlock<R>,GridMaskSet<R>)");
      return GridMaskSet<R>(A.grid(), regular_intersection(A.lattice_set(),B.lattice_set()));
    }

    template<class R>
    GridMaskSet<R>
    regular_intersection(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"regular_intersection(GridMaskSet<R>,GridMaskSet<R>)");
      if(A.block()!=B.block()) { throw IncompatibleGrids(__PRETTY_FUNCTION__); }
      return GridMaskSet<R>(A.grid(), regular_intersection(A.lattice_set(),B.lattice_set()));
    }

    template<class R>
    GridCellListSet<R>
    regular_intersection(const GridCellListSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"regular_intersection(GridCellListSet<R>,GridMaskSet<R>)");
      return GridCellListSet<R>(A.grid(), regular_intersection(A.lattice_set(),B.lattice_set()));
    }

    template<class R>
    GridCellListSet<R>
    regular_intersection(const GridMaskSet<R>& A, const GridCellListSet<R>& B)
    {
      check_same_grid(A,B,"regular_intersection(GridMaskSet<R>,GridCellListSet<R>)");
      return GridCellListSet<R>(A.grid(), regular_intersection(A.lattice_set(),B.lattice_set()));
    }


    template<class R>
    GridMaskSet<R>
    difference(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"difference(GridMaskSet<R>,GridMaskSet<R>)");
      if(A.block()!=B.block()) { throw IncompatibleGrids(__PRETTY_FUNCTION__); }
      return GridMaskSet<R>(A.grid(), difference(A.lattice_set(),B.lattice_set()));
    }

    template<class R>
    GridCellListSet<R>
    difference(const GridCellListSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"difference(GridCellListSet<R>,GridMaskSet<R>)");
      return GridCellListSet<R>(A.grid(),Combinatoric::difference(A.lattice_set(),B.lattice_set()));
    }

    template<class R>
    GridMaskSet<R>::operator ListSet<R,Rectangle>() const
    {
      check_bounded(*this);
      ListSet<R,Rectangle> result(this->dimension());
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
      if(verbosity>7) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
      GridCellListSet<R> result(g);
      check_equal_dimensions(z,g,"over_approximation(Zonotope<R>,Grid<R>)");
      if(z.empty()) {
        return result; 
      }
      Rectangle<R> bb=z.bounding_box();

      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(!disjoint(Rectangle<R>(cell),z)) {
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
      check_equal_dimensions(p,g,"over_approximation(Polytope<R>,Grid<R>)");
      if(p.empty()) {
        return result; 
      }
      Rectangle<R> bb=p.bounding_box();

      GridBlock<R> gbb=over_approximation(bb,g);
      Combinatoric::LatticeBlock block=gbb.lattice_set();

      for(Combinatoric::LatticeBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(!disjoint(Rectangle<R>(cell),p)) {
          result.adjoin(cell);
        }
      }

      return result;
    }

   
    template<class R>
    GridCellListSet<R>
    under_approximation(const Parallelotope<R>& p, const Grid<R>& g)
    {
      GridCellListSet<R> result(g);
      check_equal_dimensions(p,g,"under_approximation(Parallelotope<R>,Grid<R>)");
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
    under_approximation(const Zonotope<R>& z, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      check_equal_dimensions(z,g,"under_approximation(Zonotope<R>,Grid<R>)");
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
      check_equal_dimensions(p,g,"under_approximation(Polytope<R>,Grid<R>)");
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
      check_equal_dimensions(p,g,"under_approximation(Polyhedron<R>,Grid<R>)");
      
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
      check_equal_dimensions(gms,g,"over_approximation(GridMaskSet<R>,FiniteGrid<R>)");

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
      check_equal_dimensions(pts,g,"over_approximation(PartitionTreeSet<R>,FiniteGrid<R>)");

      for(typename PartitionTreeSet<R>::const_iterator iter=pts.begin(); iter!=pts.end(); ++iter) {
        result.adjoin(over_approximation(Rectangle<R>(*iter),g.grid()));
      }
      return result;
    }
    
    template<class R>
    GridMaskSet<R>
    over_approximation(const Set<R>& set, const FiniteGrid<R>& fg) 
    {
      GridMaskSet<R> result(fg);
      check_equal_dimensions(set,fg,"over_approximation(PartitionTreeSet<R>,FiniteGrid<R>)");
      
      const Grid<R>& g=fg.grid();
      const Combinatoric::LatticeBlock& lb=fg.bounds();
      
      for(typename Combinatoric::LatticeBlock::const_iterator iter=lb.begin(); iter!=lb.end(); ++iter) {
        GridCell<R> gc(g,*iter);
        if(!bool(set.disjoint(Rectangle<R>(gc)))) {
          result.adjoin(gc);
        }
      }
      return result;
    }
    
    // FIXME: Add GridBlock iterator
    
    template<class R>
    GridMaskSet<R>
    over_approximation(const Set<R>& set, const Grid<R>& g) 
    {
      FiniteGrid<R> fg(g,set.bounding_box());
      return over_approximation(set,fg);
    }
      
    template<class R>
    GridMaskSet<R>
    under_approximation(const ListSet<R,Rectangle>& ls, const FiniteGrid<R>& g) 
    {
      return under_approximation(GridMaskSet<R>(ls),g);
    }

    template<class R>
    GridMaskSet<R>
    under_approximation(const GridMaskSet<R>& gms, const FiniteGrid<R>& g) 
    {
      check_equal_dimensions(gms,g,"under_approximation(GridMaskSet<R>,FiniteGrid<R>)");
      
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
      return under_approximation(ListSet<R,Rectangle>(pts),g);
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
      os << "GridCellListSet<" << name<R>() << ">(\n";
      os << "  grid=" << this->grid() << "\n";
      os << "  lattice_set=" << this->lattice_set();
      os << ")" << std::endl;
      return os;
    }



    template<class R>
    std::ostream& 
    GridMaskSet<R>::write(std::ostream& os) const 
    {
      os << "GridMaskSet<" << name<R>() << ">("<< std::endl;
      os << "  grid=" << this->grid();
      os << "  bounds=" << this->bounds() << std::endl;
      os << "  mask=" << this->mask() << std::endl;
      os << ")\n";
      return os;
    }

  }
}
