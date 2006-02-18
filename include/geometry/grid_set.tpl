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

#include "../geometry/parallelopiped.h"
#include "../geometry/partition_tree_set.h"
#include "../geometry/grid_set.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    FiniteGrid<R>*
    FiniteGrid<R>::clone() const
    {
      std::cerr << "WARNING: Cloning FiniteGrid<R> causes memory leak" << std::endl;
      return new FiniteGrid<R>(*this);
    }

    template<typename R>
    FiniteGrid<R>::~FiniteGrid()
    {
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const Rectangle<R>& r, size_type n)
      : _subdivision_coordinates(r.dimension())
    {
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        R lower(r.lower_bound(i));
        R upper(r.upper_bound(i));
        R step=(upper-lower)/n;
        _subdivision_coordinates[i].push_back(lower);
        for(size_type j=1; j!=n; ++j) {
          _subdivision_coordinates[i].push_back(lower+j*step);
        }
        _subdivision_coordinates[i].push_back(upper);
      }
      create();
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const Rectangle<R>& r, SizeArray sz)
      : _subdivision_coordinates(r.dimension())
    {
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        R lower(r.lower_bound(i));
        R upper(r.upper_bound(i));
        R step=(upper-lower)/sz[i];
        _subdivision_coordinates[i].push_back(lower);
        for(size_type j=1; j!=sz[i]; ++j) {
          _subdivision_coordinates[i].push_back(lower+j*step);
        }
        _subdivision_coordinates[i].push_back(upper);
      }
      create();
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const array< std::vector<R> >& sp)
      : _subdivision_coordinates(sp)
    {
      create();
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const ListSet<R,Rectangle>& ls)
      : _subdivision_coordinates(ls.dimension())
    {
      for(typename ListSet<R,Rectangle>::const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        for(dimension_type n=0; n!=ls.dimension(); ++n) {
          _subdivision_coordinates[n].push_back(riter->lower_bound(n));
          _subdivision_coordinates[n].push_back(riter->upper_bound(n));
        }
      }
      create();
    }

    template<typename R>
    FiniteGrid<R>::FiniteGrid(const FiniteGrid<R>& g1, FiniteGrid<R>& g2)
      : _subdivision_coordinates(g1.dimension())
    {
      for(dimension_type d=0; d!=dimension(); ++d) {
        std::vector<R>& sc(_subdivision_coordinates[d]);
        const std::vector<R>& sc1(g1._subdivision_coordinates[d]);
        const std::vector<R>& sc2(g2._subdivision_coordinates[d]);
        sc.resize(sc1.size()+sc2.size());
        std::merge(sc1.begin(),sc1.end(),sc2.begin(),sc2.end(),sc.begin());
      }
      create();
    }

    template<typename R>
    void
    FiniteGrid<R>::create()
    {
       for(dimension_type i=0; i!=dimension(); ++i) {
        std::vector<R>& pos=_subdivision_coordinates[i];
        std::sort(pos.begin(),pos.end());
        typename std::vector<R>::iterator newend=std::unique(pos.begin(),pos.end());
        pos.resize(std::distance(pos.begin(),newend));
      }
      _strides=compute_strides(sizes());
    }

    template<typename R>
    array< std::vector<index_type> >
    FiniteGrid<R>::index_translation(const FiniteGrid<R>& from, const FiniteGrid<R>& to)
    {
      assert(from.dimension()==to.dimension());
      array< std::vector<index_type> > result(from.dimension());
      for(dimension_type d=0; d!=from.dimension(); ++d) {
        for(size_type n=0; n!=from.size(d); ++n) {
          index_type i=to.subdivision_index(d,from.subdivision_coordinate(d,n));
          result[d].push_back(i);
        }
      }
      return result;
    }


    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const Grid<R>& g)
    {
      return g.write(os);
    }

    template<typename R>
    std::ostream&
    FiniteGrid<R>::write(std::ostream& os) const
    {
      return os << "FiniteGrid(" << this->_subdivision_coordinates << ")";
    }

    template<typename R>
    std::ostream&
    InfiniteGrid<R>::write(std::ostream& os) const
    {
        return os << "InfiniteGrid";
    }




    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g)
      : _grid(g), _position(g.dimension())
    { }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const IndexBlock& b)
      : _grid(g), _position(b)
    {
      assert(g.dimension()==b.dimension());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid(g), _position(l,u)
    {
      assert(g.dimension()==l.size());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const Rectangle<R>& r)
      : _grid(g), _position(g.dimension())
    {
      assert(g.dimension()==r.dimension());
      for(dimension_type i=0; i!=dimension(); ++i) {
        /* TODO: Catch and rethrow exceptions */
        _position.set_lower_bound(i,g.subdivision_index(i,r.lower_bound(i)));
        _position.set_upper_bound(i,g.subdivision_index(i,r.upper_bound(i)));
      }
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const GridCell<R>& c)
      : _grid(c._grid), _position(c.dimension())
    {
      for(dimension_type i=0; i!=_position.dimension(); ++i) {
        _position.set_lower_bound(i,c.position()[i]);
        _position.set_upper_bound(i,c.position()[i]+1);
      }
    }


    template<typename R>
    GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
      : _grid(g), _position(pos)
    {
      assert(_position.size()==_grid.dimension());
    }



    template<typename R>
    GridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(dimension_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_position[i]));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_position[i]+1));
      }

      return result;
    }


    template<typename R>
    GridRectangle<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(size_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_position.lower_bound(i)));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_position.upper_bound(i)));
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
    interiors_intersect(const Rectangle<R>& A, const GridMaskSet<R>& B) {
      assert(A.dimension()==B.dimension());
      GridRectangle<R> gr=over_approximation_of_intersection(A,B.bounding_box(),B.grid());
      return interiors_intersect(gr,B);
    }
    
    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const GridRectangle<R>& A, const GridMaskSet<R>& B) {
      GridRectangleListSet<R> grls(A.grid());
      grls.push_back(A);
      GridMaskSet<R> gms(grls);
      return !regular_intersection(gms,B).empty();
    }
    

    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect (const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return !regular_intersection(A,B).empty();
    }

    /*! \brief Tests if A is a subset of the interior of B. */
    template<typename R>
    bool inner_subset (const GridMaskSet<R> &A, const GridMaskSet<R> &B);

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset (const Rectangle<R> &r, const GridMaskSet<R> &B)
    {
      assert(r.dimension() == B.dimension());
      GridMaskSet<R> A(B.grid(),B.bounds());
      A.adjoin(over_approximation(r,A.grid()));
      return subset(A,B);
    }

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset (const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return A._mask <= B._mask;
    }

    /*! \brief The union of \a A and \a B. */
    template<typename R>
    GridMaskSet<R>
    join(const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return GridMaskSet<R>(A.grid(), A.bounds(), A._mask | B._mask);
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    regular_intersection(const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return GridMaskSet<R>(A.grid(), A.bounds(), A._mask & B._mask);
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    difference(const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return GridMaskSet<R>(A.grid(), A.bounds(), A._mask - B._mask);
    }




    
    template<typename R>
    GridMaskSetConstIterator<R>::GridMaskSetConstIterator(const GridMaskSet<R>& s, 
                                                          const size_type& ind)
      : _set(s), _index(ind) 
    { 
      initialise(); 
    }

    template<typename R>
    GridMaskSetConstIterator<R>::GridMaskSetConstIterator(const GridMaskSet<R>& s,
                                                          const IndexArray& pos)
      : _set(s),
        _index(compute_index(pos,s._lower,s._strides))
    {
      initialise();
    }

    template<typename R>
    void 
    GridMaskSetConstIterator<R>::initialise() { 
      if(_index>=_set._mask.size()) {
        _index=_set._mask.size();
      }
      else if(_set._mask[_index]==false) { 
        ++(*this); 
      }
    }
    
    template<typename R>
    GridCell<R>
    GridMaskSetConstIterator<R>::operator*() const
    {
      return GridCell<R>(_set.grid(),compute_position(_index,_set._lower,_set._strides));
    }

    template<typename R>
    GridMaskSetConstIterator<R>&
    GridMaskSetConstIterator<R>::operator++()
    {
      do {
        ++_index;
      } while( _index<_set._mask.size() && _set._mask[_index]==false );
      return *this;
    }





    template<typename R>
    GridMaskSet<R>::operator ListSet<R,Rectangle>() const
    {
      //std::cerr << "GridMaskSet<R>::operator ListSet<R,Rectangle>() const" << std::endl;
      ListSet<R,Rectangle> result(dimension());
      for(typename GridMaskSet::const_iterator riter=begin(); riter!=end(); ++riter) {
        result.push_back(*riter);
      }
      return result;
    }


    template<typename R>
    GridCellListSet<R>::operator ListSet<R,Rectangle>() const
    {
      //std::cerr << "GridCellListSet<R>::operator ListSet<R,Rectangle>() const" << std::endl;
      ListSet<R,Rectangle> result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }

    template<typename R>
    GridRectangleListSet<R>::operator ListSet<R,Rectangle>() const
    {
      //std::cerr << "GridRectangleListSet<R>::operator ListSet<R,Rectangle>() const" << std::endl;
      ListSet<R,Rectangle> result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }




    template<typename R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& g)
      : _grid_ptr(&g), _bounds(g.bounds()), 
        _lower(g.lower()), _sizes(g.dimension()), _strides(g.dimension()+1),
        _mask()
    {
      //std::cerr << "GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& g)" << std::endl;
      _sizes=_bounds.sizes();
      _strides=_bounds.strides();
      _mask=BooleanArray(capacity(),false);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const IndexBlock& b)
      : _grid_ptr(&g), _bounds(b),
        _lower(g.dimension()), _sizes(g.dimension()), _strides(_sizes.size()+1),
        _mask()
    {
      _lower=_bounds.lower();
      _sizes=_bounds.sizes();
      _strides=_bounds.strides();
      _mask=BooleanArray(capacity(),false);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const IndexBlock& b, const BooleanArray& m)
      : _grid_ptr(&g), _bounds(b),
        _lower(g.dimension()), _sizes(g.dimension()), _strides(_sizes.size()+1),
        _mask(m)
    {
      _lower=_bounds.lower();
      _sizes=_bounds.sizes();
      _strides=_bounds.strides();
      assert(_mask.size()==capacity());
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
      : _grid_ptr(&gms.grid()), _bounds(gms._bounds),
        _lower(gms._lower), _sizes(gms._sizes), _strides(gms._strides),
        _mask(gms._mask)
    {  
    }
    
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()), _bounds(grid().dimension()),
        _lower(grid().dimension()), _sizes(grid().dimension()), _strides(grid().dimension()+1),
        _mask()
    {
      const FiniteGrid<R>* finite_grid_ptr=dynamic_cast<const FiniteGrid<R>*>(_grid_ptr);
      if(finite_grid_ptr) {
        _bounds=finite_grid_ptr->bounds();
      }
      else {
        _bounds=compute_cell_list_bounds(gcls._list);
      }
      _lower=_bounds.lower();
      _sizes=_bounds.sizes();
      _strides=_bounds.strides();
      _mask=BooleanArray(capacity(),false);
      compute_cell_list_mask(&_mask, _strides, _lower, gcls._list);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridRectangleListSet<R>& grls)
      : _grid_ptr(&grls.grid()), _bounds(grid().dimension()),
        _lower(grid().dimension()), _sizes(grid().dimension()), _strides(grid().dimension()+1),
        _mask()
    {
      const FiniteGrid<R>* finite_grid_ptr=dynamic_cast<const FiniteGrid<R>*>(_grid_ptr);
      if(finite_grid_ptr) {
        _bounds=finite_grid_ptr->bounds();
      }
      else {
        _bounds=compute_rectangle_list_bounds(grls._list);
      }
      _lower=_bounds.lower();
      _sizes=_bounds.sizes();
      _strides=_bounds.strides();
      _mask=BooleanArray(capacity(),false);
      compute_rectangle_list_mask(&_mask, _strides, _lower, grls._list);
    }


    template<typename R>
    GridMaskSet<R>::GridMaskSet(const ListSet<R,Rectangle>& rls) 
      : _grid_ptr(0), _bounds(rls.dimension()),
        _lower(rls.dimension()), _sizes(rls.dimension()), _strides(rls.dimension()+1),
        _mask()
    {
      //FIXME: Memory leak!    
      _grid_ptr=new FiniteGrid<R>(rls);
      _bounds=dynamic_cast<const FiniteGrid<R>*>(_grid_ptr)->bounds();
      _lower=_bounds.lower();
      _sizes=_bounds.sizes();
      _strides=_bounds.strides();
      _mask=BooleanArray(capacity(),false);

      for(typename ListSet<R,Rectangle>::const_iterator riter=rls.begin(); 
          riter!=rls.end(); ++riter) 
      {
        GridRectangle<R> r(grid(),*riter);
        adjoin(r);
      }
    }    
    
  



    template<typename R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
      : _grid_ptr(&g), _list(g.dimension())
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
      : _grid_ptr(&gms.grid()), _list(gms.dimension())
    {
      append_to_cell_list(&_list, gms._lower, gms._strides, gms._mask);
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()), _list(gcls._list)
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridRectangleListSet<R>& grls)
      : _grid_ptr(&grls.grid()), _list(grls.dimension())
    {
      append_to_cell_list(&_list, grls._list);
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const ListSet<R,Rectangle>& rls)
      : _grid_ptr(0), _list(rls.dimension())
    {
      (*this)=GridCellListSet<R>(GridRectangleListSet<R>(rls));      
    }




    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g)
      : _grid_ptr(&g), _list(g.dimension())
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const ListSet<R,Rectangle>& s)
      : _grid_ptr(new FiniteGrid<R>(s)), _list(s.dimension())
    {
      //std::cerr << "GridRectangleListSet<R>::GridRectangleListSet(const ListSet<R,Rectangle>& s)" << std::endl;
      typedef typename ListSet<R,Rectangle>::const_iterator list_set_const_iterator;
      for(list_set_const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
        this->push_back(GridRectangle<R>(grid(),*iter));
      }
    }

    /* FIXME: This constructor is only included since Boost Python doesn't find the conversion operator */
    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const PartitionTreeSet<R>& pts)
      : _grid_ptr(0),_list(pts.dimension()) 
    {
      //std::cerr << "GridRectangleListSet<R>::GridRectangleListSet(const PartitionTreeSet<R>& pts)" << std::endl;
      _grid_ptr=new FiniteGrid<R>(pts.bounding_box(),pts.subdivisions());
      for(typename PartitionTreeSet<R>::const_iterator iter=pts.begin(); iter!=pts.end(); ++iter) {
        Rectangle<R> rect(*iter);
        GridRectangle<R> grect(this->grid(),rect);
        this->push_back(grect);
      }
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const GridRectangleListSet<R>& s)
      : _grid_ptr(&s.grid()), _list(s._list)
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridCellListSet<R>& s)
      : _grid_ptr(&g), _list(g.dimension())
    {
      assert(g.dimension()==s.dimension());
      _list.resize(2*s.size());

      const FiniteGrid<R>* to=dynamic_cast<const FiniteGrid<R>*>(&g);
      const FiniteGrid<R>* from=dynamic_cast<const FiniteGrid<R>*>(&s.grid());

      array< std::vector<index_type> > tr = FiniteGrid<R>::index_translation(*from,*to);
      translate_cell_coordinates(&_list, s._list, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridRectangleListSet<R>& s)
      : _grid_ptr(&g), _list(s._list)
    {
      const FiniteGrid<R>* to=dynamic_cast<const FiniteGrid<R>*>(&g);
      const FiniteGrid<R>* from=dynamic_cast<const FiniteGrid<R>*>(&s.grid());

      array< std::vector<index_type> > tr = FiniteGrid<R>::index_translation(*from,*to);
      translate_rectangle_coordinates(&_list, s._list, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const ListSet<R,Rectangle>& ls)
      : _grid_ptr(&g), _list(ls.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator ListSet_const_iterator;

      for(ListSet_const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        _list.push_back(grid().lower_index(*riter));
        _list.push_back(grid().upper_index(*riter));
      }
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
    GridCellListSet<R>
    over_approximation(const Parallelopiped<R>& p, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(g.dimension()==p.dimension());
      if(p.empty()) {
        return result; 
      }
      Rectangle<R> bb=p.bounding_box();
      GridRectangle<R> gbb=over_approximation(bb,g);
      IndexBlock block=gbb.position();

      for(IndexBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(!disjoint(Rectangle<R>(cell),p)) {
          result.adjoin(cell);
        }
      }
      return result;
    }

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) 
    {
      //std::cerr << "GridMaskSet<R>::over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) " << std::endl;
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==ls.dimension());

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation(*iter,g));
      }
        
      return result;
    }

    template<typename R>
    GridRectangle<R>
    over_approximation_of_intersection(const Rectangle<R>& r1, 
                                       const Rectangle<R>& r2,
                                       const Grid<R>& g) 
    {
      return over_approximation(intersection(r1,r2),g);
    }

    template<typename R>
    GridCellListSet<R>
    over_approximation_of_intersection(const Parallelopiped<R>& p, 
                                       const Rectangle<R>& r,
                                       const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(p.dimension()==r.dimension());
      assert(g.dimension()==p.dimension());
      Rectangle<R> bb=intersection(p.bounding_box(),r);

      if(!bb.empty()) {
        GridRectangle<R> gbb=over_approximation(bb,g);
        IndexBlock block=gbb.position();

        for(IndexBlock::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
          GridCell<R> cell(g,*iter);
          if(!disjoint(Rectangle<R>(cell),p)) {
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
      //std::cerr << "GridMaskSet<R>::over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) " << std::endl;
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==ls.dimension());

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation_of_intersection(*iter,r,g));
      }
        
      return result;
    }



    template <typename R>
    std::ostream& operator<<(std::ostream &os,
                             const GridRectangleListSet<R>& set)
    {
      os << "GridRectangleListSet<" << name<R>() << ">(\n  rectangle_list: [\n    ";
      for (size_type i=0; i!=set.size(); i++) {
        if(i!=0) {
          os << ",\n    ";
        }
        os << Rectangle<R>(set[i]);
      }
      os << "\n  ]\n";
      os << "  grid: " << set.grid() << "\n";
      os << "  integer_rectangle_list: [ ";
      os.flush();
      for(size_type i=0; i!=set.size(); ++i) {
        if(i!=0) { os << ", "; }
        os << set[i].position();
      }
      os << " ]\n";
      os << ")" << std::endl;
      return os;
    }



    template <typename R>
    std::ostream& operator<<(std::ostream &os,
                             const GridCellListSet<R>& set)
    {
      os << "GridCellListSet<" << name<R>() << ">(\n";
      /*
      os << "  rectangle_list: [\n    ";
      if (!set.empty() ) {
        os << set[0];
      }
      for (size_t i=1; i<set.size(); ++i) {
        os << ",\n    " << Rectangle<R>(set[i]);
      }
      os << "\n  ]\n";
      */
      os << "  grid: " << set.grid() << "\n";
      os << "  integer_cell_list: [ ";
      os.flush();
      for(size_type i=0; i!=set.size(); ++i) {
        if(i!=0) { os << ", "; }
        os <<  set[i].position();
      }
      os << " ]\n";
      os << ")" << std::endl;
      return os;
    }

    template <typename R>
    std::ostream& operator<<(std::ostream &os,
                             const GridMaskSet<R>& set)
    {
      os << "GridMaskSet<" << name<R>() << ">(\n";
      os << "  grid: " << set.grid() << "\n";
      os << "  mask: " << set.mask() << "\n";
      os << ")\n";
      return os;
    }

  }
}
