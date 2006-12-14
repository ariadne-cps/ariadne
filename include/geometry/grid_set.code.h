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

#include "../geometry/ppl_polyhedron.h"
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
    

    
    // GridBlock --------------------------------------------------------------
    
    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g)
      : _grid(g), _lattice_set(g.dimension())
    { 
      _lattice_set.set_lower_bound(0,1);
      //_lattice_set.set_lower_bound(0,0);
      _lattice_set.set_upper_bound(0,0);
    }

    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const Combinatoric::LatticeBlock& b)
      : _grid(g), _lattice_set(b)
    {
      check_equal_dimensions(g,b,"GridBlock<R>::GridBlock(Grid<R>,LatticeBlock)");
    }

    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid(g), _lattice_set(l,u)
    {
      check_dimension(g,l.size(),__PRETTY_FUNCTION__);
    }

    template<class R>
    GridBlock<R>::GridBlock(const Grid<R>& g, const Rectangle<R>& r)
      : _grid(g), _lattice_set(g.dimension())
    {
      check_equal_dimensions(g,r,"GridBlock<R>::GridBlock(Grid<R>,Rectangle<R>)");
      for(dimension_type i=0; i!=dimension(); ++i) {
        /* TODO: Catch and rethrow exceptions */
        _lattice_set.set_lower_bound(i,g.subdivision_index(i,r.lower_bound(i)));
        _lattice_set.set_upper_bound(i,g.subdivision_index(i,r.upper_bound(i)));
      }
    }

    template<class R>
    GridBlock<R>::GridBlock(const GridCell<R>& c)
      : _grid(c._grid), _lattice_set(c.lattice_set())
    {
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
    subset(const GridBlockListSet<R>& A, const GridMaskSet<R>& B)
    {
      check_same_grid(A,B,"subset(GridBlockListSet<R>,GridMaskSet<R>)");
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

    template<class R>
    GridCellListSet<R>::GridCellListSet(const GridBlockListSet<R>& grls)
      : _grid_ptr(&grls.grid()), _lattice_set(grls.dimension())
    {
      this->_lattice_set.adjoin(grls._lattice_set); 
    }

    template<class R>
    GridCellListSet<R>::GridCellListSet(const ListSet<R,Rectangle>& rls)
      : _grid_ptr(0), _lattice_set(rls.dimension())
    {
      (*this)=GridCellListSet<R>(GridBlockListSet<R>(rls));      
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

    // GridBlockListSet ------------------------------------------------------------

    template<class R>
    GridBlockListSet<R>::GridBlockListSet(const Grid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
    }

    template<class R>
    GridBlockListSet<R>::GridBlockListSet(const Grid<R>& g, 
                                        const Combinatoric::LatticeBlockListSet& lbls)
      : _grid_ptr(&g), _lattice_set(lbls)
    {
      check_equal_dimensions(g,lbls,"GridBlockListSet<R>::GridBlockListSet(Grid<R>,LatticeBlockListSet)");
    }

    // FIXME: Memory leak
    template<class R>
    GridBlockListSet<R>::GridBlockListSet(const ListSet<R,Rectangle>& s)
      : _grid_ptr(new IrregularGrid<R>(s)), _lattice_set(s.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator list_set_const_iterator;
      for(list_set_const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
        this->adjoin(GridBlock<R>(grid(),*iter));
      }
    }

    /* FIXME: This constructor is only included since Boost Python doesn't find the conversion operator */
    // FIXME: Memory leak due to construction of new grid
    template<class R>
    GridBlockListSet<R>::GridBlockListSet(const PartitionTreeSet<R>& pts)
      : _grid_ptr(0),_lattice_set(pts.dimension()) 
    {
      SizeArray sizes=pts.depths();
      for(dimension_type i=0; i!=pts.dimension(); ++i) {
        sizes[i]=exp2(sizes[i]);
      }
      _grid_ptr=new IrregularGrid<R>(pts.bounding_box(),sizes);
      for(typename PartitionTreeSet<R>::const_iterator iter=pts.begin(); iter!=pts.end(); ++iter) {
        Rectangle<R> rect(*iter);
        GridBlock<R> grect(this->grid(),rect);
        this->adjoin(grect);
      }
    }

    template<class R>
    GridBlockListSet<R>::GridBlockListSet(const GridBlockListSet<R>& s)
      : _grid_ptr(&s.grid()), _lattice_set(s._lattice_set)
    {
    }

    
    template<class R>
    GridBlockListSet<R>::GridBlockListSet(const Grid<R>& g, const ListSet<R,Rectangle>& ls)
      : _grid_ptr(&g), _lattice_set(ls.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator ListSet_const_iterator;

      for(ListSet_const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        _lattice_set.adjoin(Combinatoric::LatticeBlock(grid().lower_index(*riter),grid().upper_index(*riter)));
      }
    }

    template<class R>
    GridBlockListSet<R>::operator ListSet<R,Rectangle>() const
    {
      ListSet<R,Rectangle> result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }

    template<class R>
    void
    GridBlockListSet<R>::clear()
    {
      this->_lattice_set.clear();
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
        assert(subset(b,irregular_grid_ptr->block()));
      }
    }

    template<class R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeMaskSet& ms)
      : _grid_ptr(&g), _lattice_set(ms)
    {
      const IrregularGrid<R>* irregular_grid_ptr=dynamic_cast<const IrregularGrid<R>*>(_grid_ptr);
      if(irregular_grid_ptr) {
        assert(subset(ms.block(),irregular_grid_ptr->block()));
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
    
    template<class R>
    GridMaskSet<R>::GridMaskSet(const GridBlockListSet<R>& grls)
      : _grid_ptr(&grls.grid()),
        _lattice_set(grls.lattice_set())
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
      check_equal_dimensions(p,g,"over_approximation(Zonotope<R>,Grid<R>)");
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

    
    
    template<class R, template<class> class BS>
    GridMaskSet<R>
    over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      
      check_equal_dimensions(ls,g,"over_approximation(ListSet<R,BS>,FiniteGrid<R>)");

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation(*iter,g.grid()));
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
      return under_approximation(GridMaskSet<R>(GridBlockListSet<R>(pts)),g);
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
    operator<<(std::ostream& os, const GridCell<R>& c) {
      return os << Rectangle<R>(c);
    }

    
    template<class R>
    std::ostream&
    operator<<(std::ostream& os, const GridBlock<R>& r) {
      return os << Rectangle<R>(r);
    }

    
    template<class R>
    std::ostream& operator<<(std::ostream& os,
                             const GridCellListSet<R>& set)
    {
      os << "GridCellListSet<" << name<R>() << ">(\n";
      /*
      os << "  rectangle_lattice_set: [\n    ";
      if (!set.empty() ) {
        os << set[0];
      }
      for (size_t i=1; i<set.size(); ++i) {
        os << ",\n    " << Rectangle<R>(set[i]);
      }
      os << "\n  ]\n";
      */
      os << "  grid=" << set.grid() << "\n";
      os << "  integer_cell_lattice_set=[ ";
      os.flush();
      for(size_type i=0; i!=set.size(); ++i) {
        if(i!=0) { os << ", "; }
        os <<  set[i].lattice_set();
      }
      os << " ]\n";
      os << ")" << std::endl;
      return os;
    }

    template<class R>
    std::ostream& operator<<(std::ostream& os,
                             const GridBlockListSet<R>& set)
    {
      os << "GridBlockListSet<" << name<R>() << ">(\n  rectangle_lattice_set=[\n    ";
      for (size_type i=0; i!=set.size(); i++) {
        if(i!=0) {
          os << ",\n    ";
        }
        os << Rectangle<R>(set[i]);
      }
      os << "\n  ]\n";
      os << "  grid=" << set.grid() << "\n";
      os << "  integer_rectangle_lattice_set=[ ";
      os.flush();
      for(size_type i=0; i!=set.size(); ++i) {
        if(i!=0) { os << ", "; }
        os << set[i].lattice_set();
      }
      os << " ]\n";
      os << ")" << std::endl;
      return os;
    }


    template<class R>
    std::ostream& operator<<(std::ostream& os,
                             const GridMaskSet<R>& set)
    {
      os << "GridMaskSet<" << name<R>() << ">("<< std::endl;
      os << "  grid=" << set.grid();
      os << "  bounds=" << set.bounds() << std::endl;
      os << "  mask=" << set.mask() << std::endl;
      os << ")\n";
      return os;
    }

  }
}
