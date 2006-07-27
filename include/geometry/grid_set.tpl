/***************************************************************************
 *            grid_set.tpl
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

#include "../utility/stlio.h"

#include "../base/array_operations.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/list_set.h"
#include "../geometry/partition_tree_set.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g)
      : _grid(g), _lattice_set(g.dimension())
    { 
      _lattice_set.set_lower_bound(0,1);
      _lattice_set.set_upper_bound(0,0);
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const LatticeRectangle& b)
      : _grid(g), _lattice_set(b)
    {
      assert(g.dimension()==b.dimension());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid(g), _lattice_set(l,u)
    {
      assert(g.dimension()==l.size());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const Rectangle<R>& r)
      : _grid(g), _lattice_set(g.dimension())
    {
      assert(g.dimension()==r.dimension());
      for(dimension_type i=0; i!=dimension(); ++i) {
        /* TODO: Catch and rethrow exceptions */
        _lattice_set.set_lower_bound(i,g.subdivision_index(i,r.lower_bound(i)));
        _lattice_set.set_upper_bound(i,g.subdivision_index(i,r.upper_bound(i)));
      }
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const GridCell<R>& c)
      : _grid(c._grid), _lattice_set(c.lattice_set())
    {
    }


    template<typename R>
    GridCell<R>::GridCell(const Grid<R>& g, const LatticeCell& pos)
      : _grid(g), _lattice_set(pos)
    {
      assert(_lattice_set.dimension()==_grid.dimension());
    }

    template<typename R>
    GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
      : _grid(g), _lattice_set(pos)
    {
      assert(_lattice_set.dimension()==_grid.dimension());
    }



    template<typename R>
    GridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(dimension_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_lattice_set.lower_bound(i)));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_lattice_set.upper_bound(i)));
      }

      return result;
    }


    template<typename R>
    GridRectangle<R>::operator Rectangle<R>() const {
      if(this->empty()) {
        return Rectangle<R>();
      }
     

      Rectangle<R> result(dimension());
      for(size_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_lattice_set.lower_bound(i)));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_lattice_set.upper_bound(i)));
      }

      return result;
    }


    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const GridRectangle<R>& r) {
      return os << Rectangle<R>(r);
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const GridCell<R>& c) {
      return os << Rectangle<R>(c);
    }



    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const GridRectangle<R>& A, const GridRectangle<R>& B) {
      assert(A.grid() == B.grid());
      return interiors_intersect(A.lattice_set(),B.lattice_set());
    }

    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const GridRectangle<R>& A, const GridMaskSet<R>& B) {
      assert(A.grid() == B.grid());
      return interiors_intersect(A._lattice_set,B._lattice_set);
    }

    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid());
      return interiors_intersect(A._lattice_set,B._lattice_set);
    }



    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const GridRectangle<R>& A, const GridRectangle<R>& B)
    {
      assert(A.dimension() == B.dimension());
      if(A.grid()==B.grid()) {
        return subset(A.lattice_set(),B.lattice_set());
      }
      return subset(Rectangle<R>(A),Rectangle<R>(B));
    }

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const GridRectangle<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid());
      return subset(A.lattice_set(),B.lattice_set()); 
    }

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid());
      return subset(A.lattice_set(),B.lattice_set());
    }



    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const Rectangle<R>& A, const GridMaskSet<R>& B) {
      assert(A.dimension()==B.dimension());
      GridRectangle<R> gr=over_approximation_of_intersection(A,Rectangle<R>(B.bounding_box()),B.grid());
      return interiors_intersect(gr,B);
    }
    
    /*! \brief Tests if A is a subset of the interior of B. */
    template<typename R>
    bool
    inner_subset(const Rectangle<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.dimension() == B.dimension());
      return subset(outer_approximation(A,B.grid()),B);
    }

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const Rectangle<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.dimension() == B.dimension());
      return subset(over_approximation(A,B.grid()),B);
    }



    /*! \brief The union of \a A and \a B. */
    template<typename R>
    GridMaskSet<R>
    join(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid() && A.bounds()==B.bounds());
      return GridMaskSet<R>(A.grid(), A.bounds(), A.mask() | B.mask());
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    regular_intersection(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid() && A.bounds()==B.bounds());
      return GridMaskSet<R>(A.grid(), A.bounds(), A.mask() & B.mask());
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    difference(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid() && A.bounds()==B.bounds());
      return GridMaskSet<R>(A.grid(), A.bounds(), A.mask() - B.mask());
    }

    template<typename R>
    GridMaskSet<R>::operator ListSet<R,Rectangle>() const
    {
      ListSet<R,Rectangle> result(this->dimension());
      for(typename GridMaskSet::const_iterator riter=begin(); riter!=end(); ++riter) {
        Rectangle<R> r(*riter);
        result.push_back(r);
      }
      return result;
    }


    template<typename R>
    GridCellListSet<R>::operator ListSet<R,Rectangle>() const
    {
      ListSet<R,Rectangle> result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }

    template<typename R>
    GridRectangleListSet<R>::operator ListSet<R,Rectangle>() const
    {
      ListSet<R,Rectangle> result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }



  
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& g)
      : _grid_ptr(&g.grid()), _lattice_set(g.bounds()) 
    { 
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeRectangle& b, const BooleanArray& m)
      : _grid_ptr(&g), _lattice_set(b,m)
    {
      const IrregularGrid<R>* irregular_grid_ptr=dynamic_cast<const IrregularGrid<R>*>(_grid_ptr);
      if(irregular_grid_ptr) {
        assert(subset(b,irregular_grid_ptr->bounds()));
      }
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeMaskSet& ms)
      : _grid_ptr(&g), _lattice_set(ms)
    {
      const IrregularGrid<R>* irregular_grid_ptr=dynamic_cast<const IrregularGrid<R>*>(_grid_ptr);
      if(irregular_grid_ptr) {
        assert(subset(ms.bounds(),irregular_grid_ptr->bounds()));
      }
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg, const BooleanArray& m)
      : _grid_ptr(&fg.grid()), _lattice_set(fg.bounds(),m)
    {
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
      : _grid_ptr(&gms.grid()), _lattice_set(gms._lattice_set)
    {
    }
    
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()),
        _lattice_set(gcls.dimension())
    {
      _lattice_set=LatticeMaskSet(gcls.lattice_set().bounds());
      _lattice_set.adjoin(gcls.lattice_set());
    }
    
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridRectangleListSet<R>& grls)
      : _grid_ptr(&grls.grid()),
        _lattice_set(grls.lattice_set().bounds())
    {
      _lattice_set.adjoin(grls.lattice_set());
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const ListSet<R,Rectangle>& rls) 
      : _grid_ptr(new IrregularGrid<R>(rls)), 
        _lattice_set(dynamic_cast<const IrregularGrid<R>*>(_grid_ptr)->bounds())
    {
      //FIXME: Memory leak!    

      for(typename ListSet<R,Rectangle>::const_iterator riter=rls.begin(); 
          riter!=rls.end(); ++riter) 
      {
        GridRectangle<R> r(grid(),*riter);
        adjoin(r);
      }
    }    
    
  



    template<typename R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
      : _grid_ptr(&gms.grid()), _lattice_set(gms.dimension())
    {
      this->_lattice_set.adjoin(gms._lattice_set);
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()), _lattice_set(gcls._lattice_set)
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridRectangleListSet<R>& grls)
      : _grid_ptr(&grls.grid()), _lattice_set(grls.dimension())
    {
      this->_lattice_set.adjoin(grls._lattice_set); 
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const ListSet<R,Rectangle>& rls)
      : _grid_ptr(0), _lattice_set(rls.dimension())
    {
      (*this)=GridCellListSet<R>(GridRectangleListSet<R>(rls));      
    }




    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const ListSet<R,Rectangle>& s)
      : _grid_ptr(new IrregularGrid<R>(s)), _lattice_set(s.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator list_set_const_iterator;
      for(list_set_const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
        this->adjoin(GridRectangle<R>(grid(),*iter));
      }
    }

    /* FIXME: This constructor is only included since Boost Python doesn't find the conversion operator */
    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const PartitionTreeSet<R>& pts)
      : _grid_ptr(0),_lattice_set(pts.dimension()) 
    {
      SizeArray sizes=pts.depths();
      for(dimension_type i=0; i!=pts.dimension(); ++i) {
        sizes[i]=pow(2,sizes[i]);
      }
      _grid_ptr=new IrregularGrid<R>(pts.bounding_box(),sizes);
      for(typename PartitionTreeSet<R>::const_iterator iter=pts.begin(); iter!=pts.end(); ++iter) {
        Rectangle<R> rect(*iter);
        GridRectangle<R> grect(this->grid(),rect);
        this->adjoin(grect);
      }
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const GridRectangleListSet<R>& s)
      : _grid_ptr(&s.grid()), _lattice_set(s._lattice_set)
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridCellListSet<R>& s)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
      //TODO: re-write this function
      assert(false);
      assert(g.dimension()==s.dimension());

      const IrregularGrid<R>* to=dynamic_cast<const IrregularGrid<R>*>(&g);
      const IrregularGrid<R>* from=dynamic_cast<const IrregularGrid<R>*>(&s.grid());

      array< std::vector<index_type> > tr = IrregularGrid<R>::index_translation(*from,*to);
      //translate_cell_coordinates(&_lattice_set, s._lattice_set, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridRectangleListSet<R>& s)
      : _grid_ptr(&g), _lattice_set(s._lattice_set)
    {
      //TODO: re-write this function
      assert(false);

      const IrregularGrid<R>* to=dynamic_cast<const IrregularGrid<R>*>(&g);
      const IrregularGrid<R>* from=dynamic_cast<const IrregularGrid<R>*>(&s.grid());

      array< std::vector<index_type> > tr = IrregularGrid<R>::index_translation(*from,*to);
      //translate_rectangle_coordinates(&_lattice_set, s._lattice_set, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const ListSet<R,Rectangle>& ls)
      : _grid_ptr(&g), _lattice_set(ls.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator ListSet_const_iterator;

      for(ListSet_const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        _lattice_set.adjoin(LatticeRectangle(grid().lower_index(*riter),grid().upper_index(*riter)));
      }
    }

    template<typename R>
    GridRectangle<R>
    outer_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      if(r.empty()) {
        return GridRectangle<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_upper_index(i,r.lower_bound(i))-1;
        upper[i]=g.subdivision_lower_index(i,r.upper_bound(i))+1;
      }
      return GridRectangle<R>(g,lower,upper);
    }
    
    template<typename R>
    GridRectangle<R>
    over_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      if(r.empty()) {
        return GridRectangle<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,r.lower_bound(i));
        upper[i]=g.subdivision_upper_index(i,r.upper_bound(i));
      }

      return GridRectangle<R>(g,lower,upper);
    }
    
    template<typename R>
    GridRectangle<R>
    under_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      
      if(r.empty()) {
        return GridRectangle<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,r.lower_bound(i))+1;
        upper[i]=g.subdivision_upper_index(i,r.upper_bound(i))-1;
      }

      return GridRectangle<R>(g,lower,upper);
    }
    
    template<typename R>
    GridCellListSet<R>
    over_approximation(const Zonotope<R>& p, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(g.dimension()==p.dimension());
      if(p.empty()) {
        return result; 
      }
      Rectangle<R> bb=p.bounding_box();

      GridRectangle<R> gbb=over_approximation(bb,g);
      LatticeRectangle block=gbb.lattice_set();

      for(LatticeRectangle::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(!disjoint(Rectangle<R>(cell),p)) {
          result.adjoin(cell);
        }
      }

      return result;
    }

    template<typename R>
    inline
    GridCellListSet<R>
    over_approximation(const Parallelotope<R>& p, const Grid<R>& g) {
       return over_approximation(Zonotope<R>(p), g);
    }
    
    template<typename R>
    GridCellListSet<R>
    under_approximation(const Zonotope<R>& z, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(g.dimension()==z.dimension());
      if(z.empty()) {
        return result; 
      }
      
      Rectangle<R> bb=z.bounding_box();
      GridRectangle<R> gbb=over_approximation(bb,g);
      LatticeRectangle block=gbb.lattice_set();

      for(LatticeRectangle::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(z.contains(Rectangle<R>(cell))) {
          result.adjoin(cell);
        }
      }
      return result;
    }

    template<typename R>
    inline
    GridCellListSet<R>
    under_approximation(const Parallelotope<R>& p, const Grid<R>& g) {
       return under_approximation(Zonotope<R>(p), g);
    }
   
    template<typename R>
    GridRectangle<R>
    over_approximation(const Rectangle<R>& r, const FiniteGrid<R>& g) 
    {
      Rectangle<R> rect=regular_intersection(r,g.bounding_box());
     
      if(rect.empty()) {
        return GridRectangle<R>(g.grid());
      }
      IndexArray lower(rect.dimension());
      IndexArray upper(rect.dimension());
      
      for(size_type i=0; i!=rect.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,rect.lower_bound(i));
        upper[i]=g.subdivision_upper_index(i,rect.upper_bound(i));
      }

      return GridRectangle<R>(g.grid(),lower,upper);
    }
    
    template<typename R>
    GridRectangle<R>
    under_approximation(const Rectangle<R>& r, const FiniteGrid<R>& g) 
    {
      Rectangle<R> rect=regular_intersection(r,g.bounding_box());
     
      if(rect.empty()) {
        return GridRectangle<R>(g.grid());
      }
      IndexArray lower(rect.dimension());
      IndexArray upper(rect.dimension());
      
      for(size_type i=0; i!=rect.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,rect.lower_bound(i))+1;
        upper[i]=g.subdivision_upper_index(i,rect.upper_bound(i))-1;
      }
       
      return GridRectangle<R>(g.grid(),lower,upper);
    }
   
    template<typename R>
    GridCellListSet<R>
    over_approximation(const Zonotope<R>& p, const FiniteGrid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(g.dimension()==p.dimension());
      if(p.empty()) {
        return result; 
      }
      
      Rectangle<R> bb=regular_intersection(Polyhedron<R>(p),
                              Polyhedron<R>(g.bounding_box())).bounding_box();
      
      if (bb.empty())
      	return result;

      GridRectangle<R> gbb=over_approximation(bb,g);
      LatticeRectangle block=gbb.lattice_set();

      for(LatticeRectangle::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {

        GridCell<R> cell(g.grid(),*iter);
        if(!disjoint(Rectangle<R>(cell),p)) {
          result.adjoin(cell);
        }
      }

      return result;
    }

    template<typename R>
    inline
    GridCellListSet<R>
    over_approximation(const Parallelotope<R>& p, const FiniteGrid<R>& g) {
       return over_approximation(Zonotope<R>(p), g);
    }
    
    template<typename R>
    GridCellListSet<R>
    under_approximation(const Zonotope<R>& p, const FiniteGrid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(g.dimension()==p.dimension());
      if(p.empty()) {
        return result; 
      }
      Rectangle<R> bb=regular_intersection(Polyhedron<R>(p),
      				Polyhedron<R>(g.bounding_box())).bounding_box();

      if (bb.empty()) {
      	return result;
      }
      
      GridRectangle<R> gbb=over_approximation(bb,g);
      LatticeRectangle block=gbb.lattice_set();

      for(LatticeRectangle::const_iterator iter=block.begin(); 
            iter!=block.end(); ++iter) {

        GridCell<R> cell(g.grid(),*iter);
        if(p.contains(Rectangle<R>(cell))) {
          result.adjoin(cell);
        }
      }

      return result;
    }

    template<typename R>
    inline
    GridCellListSet<R>
    under_approximation(const Parallelotope<R>& p, const FiniteGrid<R>& g) {
       return under_approximation(Zonotope<R>(p), g);
    }
    
    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==ls.dimension());

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation(*iter,g));
      }
        
      return result;
    }
     
    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    under_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==ls.dimension());
      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(under_approximation(*iter,g));
      }
        
      return result;
    }

    template<typename R>
    GridMaskSet<R>
    over_approximation(const GridMaskSet<R>& gm, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==gm.dimension());

      for(size_t i=0; i< gm.size(); i++) 
        result.adjoin(over_approximation(Rectangle<R>(gm[i]),g));
        
      return result;
    }
    
    template<typename R>
    GridMaskSet<R>
    under_approximation(const GridMaskSet<R>& gm, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==gm.dimension());
      for(size_t i=0; i< gm.size(); i++) {
        result.adjoin(under_approximation(Rectangle<R>(gm[i]),g));
      }
        
      return result;
    }

/* NEW CODE */ 
 
    template<typename R>
    inline
    GridRectangle<R>
    over_approximation(const Rectangle<R>& r, const GridMaskSet<R>& gms) 
    {
      return over_approximation(r,gms.grid());
    }
    
    template<typename R>
    inline
    GridRectangle<R>
    under_approximation(const Rectangle<R>& r, const GridMaskSet<R>& gms) 
    {
      return under_approximation(r,gms.grid());
    }
    
    template<typename R, template < typename R > class BS >
    GridCellListSet<R>
    over_approximation_not_intersecting(const BS<R>& bs, 
                                        const GridMaskSet<R>& gms) 
    {
      GridCellListSet<R> result(gms.grid());
      assert(gms.dimension()==bs.dimension());
      if(bs.empty()) {
        return result; 
      }

      GridMaskSet<R> new_gms(gms);
      Rectangle<R> bb_gms=Rectangle<R>(new_gms.bounding_box());
      bb_gms=regular_intersection(bs.bounding_box(),bb_gms);
      
      new_gms.adjoin(over_approximation(bb_gms,gms.grid()));
      
      new_gms=difference(new_gms,gms);

      size_type new_gms_elem=new_gms.size();      

      for (size_type i=0; i< new_gms_elem; i++) {
        GridCell<R> cell(new_gms[i]);
        if(!disjoint(Rectangle<R>(cell),bs)) {
          result.adjoin(cell);
        }
      }

      return result;
    }

    template<typename R, template < typename R > class BS >
    GridCellListSet<R>
    under_approximation_not_intersecting(const BS<R>& bs, 
                                         const GridMaskSet<R>& gms) 
    {
      GridCellListSet<R> result(gms.grid());
      assert(gms.dimension()==bs.dimension());
      if(bs.empty()) {
        return result; 
      }
      
      GridMaskSet<R> new_gms(gms);
      Rectangle<R> bb_gms=Rectangle<R>(new_gms.bounding_box());
      bb_gms=regular_intersection(bs.bounding_box(),bb_gms);
      
      new_gms.adjoin(over_approximation(bb_gms,gms.grid()));

      new_gms=difference(new_gms,gms);

      size_type new_gms_elem=new_gms.size();      

      for (size_type i=0; i< new_gms_elem; i++) {
        GridCell<R> cell(new_gms[i]);
        if(bs.contains(Rectangle<R>(cell))) {
          result.adjoin(cell);
        }
      }
      
      return result;
    }

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    join_over_approximation(const GridMaskSet<R>& gms, const BS<R>& bs) 
    {
      GridMaskSet<R> result(gms);
      
      assert(gms.dimension()==bs.dimension());

      result.adjoin(over_approximation_not_intersecting(bs,gms));
        
      return result;
    }

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    join_under_approximation(const GridMaskSet<R>& gms, const BS<R>& bs) 
    {
      GridMaskSet<R> result(gms);
      
      assert(gms.dimension()==bs.dimension());

      result.adjoin(under_approximation_not_intersecting(bs,gms));
        
      return result;
    }

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    join_over_approximation(const GridMaskSet<R>& gms, const ListSet<R,BS>& ls) 
    {
      GridMaskSet<R> result(gms);
      
      assert(gms.dimension()==ls.dimension());

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); 
                                              iter!=ls.end(); ++iter) {

       result.adjoin(over_approximation_not_intersecting(*iter,gms));
      }
        
      return result;
    }
     
    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    join_under_approximation(const GridMaskSet<R>& gms, 
                             const ListSet<R,BS>& ls) 
    {
      GridMaskSet<R> result(gms);
      
      assert(gms.dimension()==ls.dimension());
      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); 
                                               iter!=ls.end(); ++iter) {
       result.adjoin(under_approximation_not_intersecting(*iter,gms));
      }
        
      return result;
    }
 
/* END NEW CODE */ 
    
    template<typename R>
    GridRectangle<R>
    over_approximation_of_intersection(const Rectangle<R>& r1, 
                                       const Rectangle<R>& r2,
                                       const Grid<R>& g) 
    {
      return over_approximation(intersection(r1,Rectangle<R>(r2)),g);
    }

    template<typename R>
    GridCellListSet<R>
    over_approximation_of_intersection(const Parallelotope<R>& p, 
                                       const Rectangle<R>& r,
                                       const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      GridRectangle<R> gbb=over_approximation_of_intersection(p.bounding_box(),r,g);
      if(!gbb.empty()) {
        LatticeRectangle block=gbb.lattice_set();
        for(LatticeRectangle::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
          GridCell<R> cell(g,*iter);
          if(!disjoint(Rectangle<R>(cell),p)) {
            result.adjoin(cell);
          }
        }
      }
      return result;
    }

    template<typename R>
    GridCellListSet<R>
    over_approximation_of_intersection(const Zonotope<R>& z, 
                                       const Rectangle<R>& r,
                                       const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      GridRectangle<R> gbb=over_approximation_of_intersection(z.bounding_box(),r,g);

      if(!gbb.empty()) {
        LatticeRectangle block=gbb.lattice_set();
        for(LatticeRectangle::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
          GridCell<R> cell(g,*iter);
          if(!disjoint(Rectangle<R>(cell),z)) {
            result.adjoin(cell);
          }
        }
      }
      return result;
    }

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation_of_intersection(const ListSet<R,BS>& ls, const Rectangle<R>& r, const FiniteGrid<R>& g) 
    {
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==ls.dimension());

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation_of_intersection(*iter,r,g.grid()));
      }
        
      return result;
    }

  

    template <typename R>
    std::ostream& operator<<(std::ostream& os,
                             const GridRectangleListSet<R>& set)
    {
      os << "GridRectangleListSet<" << name<R>() << ">(\n  rectangle_lattice_set=[\n    ";
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



    template <typename R>
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

    template <typename R>
    std::ostream& operator<<(std::ostream& os,
                             const GridMaskSet<R>& set)
    {
      os << "GridMaskSet<" << name<R>() << ">(\n";
      os << "  grid=" << set.grid() << "\n";
      os << "  bounds=" << set.bounds() << "\n";
      os << "  mask=" << set.mask() << "\n";
      os << ")\n";
      return os;
    }

  }
}
