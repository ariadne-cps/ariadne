/***************************************************************************
 *            grid_denotable_set.h
 *
 *  10 January 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file grid_denotable_set.h
 *  \brief Denotable sets.
 */

#ifndef _GRID_DENOTABLE_SET_H
#define _GRID_DENOTABLE_SET_H

#include <iostream>

#include "array.h"
#include "interval.h"
#include "binary_word.h"
#include "rectangle.h"
#include "state.h"

#include "grid.h"
#include "grid_rectangle.h"

#include "list_set.h"

/** \internal
 * Non-templated functions needed to implement grid denotable sets.
 *
 * Grid:
 * Insert a cell into a mask array
 *    insert(BooleanArray& mask, const IntegerArray& index, const IntegerArray& strides);
 * Insert a list of cells into a mask array
 *   insert(BooleanArray& mask, const IntegerBlockArray& indices, const IntegerArray& strides);
 * Insert a rectangle into a mask array.
 *   insert_rectangle(BooleanArray& mask, const IntegerArray& lower, const IntegerArray& upper, const IntegerArray& strides);
 * Insert a list of rectangles into a mask array.
 *   insert_rectangle(BooleanArray& mask, const IntegerBlockArray& indices, const IntegerArray& strides);
 * and, or and subtraction of mask arrays
 *   &=(BooleanArray& mask1, const BooleanArray& mask2)
 *   |=(BooleanArray& mask1, const BooleanArray& mask2)
 *   inplace_and_not(BooleanArray& mask1, const BooleanArray& mask2)
 * Extract list of cells from a mask array
 * Extract a list of rectangles from a mask array
 * Compute bounds from a list of cells.
 * Compute bounds from a list of rectangles.
 */


namespace Ariadne {
  namespace Geometry {
    template<typename R, template<typename> class BS> class ListSet;
    template<typename R> class Rectangle;
    template<typename R> class State;
    template<typename R> class Grid;

    template<typename R> class GridMaskSet;
    template<typename R> class GridRectangleListSet;
    template<typename R> class GridCellListSet;

    template<typename R> GridMaskSet<R> join(const GridMaskSet<R>&, const GridMaskSet<R>&);

    template<typename R> std::ostream& operator<<(std::ostream&, const GridRectangleListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridCellListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridMaskSet<R>&);

    template<class R>
    class GridMaskSet_const_iterator {
      typedef GridMaskSet_const_iterator Self;
     public:
      GridMaskSet_const_iterator(const GridMaskSet<R>& s, const IndexArray& pos);
      GridMaskSet_const_iterator(const GridMaskSet<R>& s, const size_type& ind)
        : _set(s), _index(ind) { }
      GridCell<R> operator*() const;
      Self& operator++();
      bool operator==(const Self& other) const { return _index==other._index && (&_set)==(&(other._set)); }
      bool operator!=(const Self& other) const { return !(*this==other); }
    private:
     public:
      const GridMaskSet<R>& _set;
      size_type _index;
    };

    /*! \brief A denotable set on a finite grid, defined using a mask. */
    template<typename R>
    class GridMaskSet {
      friend class GridMaskSet_const_iterator<R>;
      friend class GridCellListSet<R>;
     public:
      typedef R real_type;
      typedef State<R> state_type;
      typedef GridMaskSet_const_iterator<R> const_iterator;

      /*!\brief Construct an empty set from a grid, and a list of lower and upper bounds. */
      GridMaskSet(const Grid<R>& g, IndexArray l, IndexArray u);

      /*!\brief Construct an empty set from a grid, a list of lower and upper bounds, and a mask */
      GridMaskSet(const Grid<R>& g, IndexArray l, IndexArray u, const BooleanArray& m);

      /*!\brief Copy constructor. */
      GridMaskSet(const GridMaskSet<R>& cls);

      /*!\brief Construct from a GridCellListSet. */
      GridMaskSet(const GridCellListSet<R>& cls);

      /*!\brief Construct from a GridRectangleListSet. */
      GridMaskSet(const GridRectangleListSet<R>& rls);
      
      /*!\brief Convert to a ListSet of Rectangles. */
      operator ListSet<R,Rectangle>() const;

      /*! \brief Equality operator. */
      bool operator==(const GridMaskSet<R>& gms) {
        if(this->_grid==gms._grid && this->_lower==gms._lower && this->_upper==gms._upper) {
          return this->_mask==gms._mask;
        }
        throw(std::domain_error("Can only compare GridMaskSets on the same grid"));
      }

      /*! \brief Inequality operator. */
      bool operator!=(const GridMaskSet<R>& gms) { return !(*this==gms); }

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return _grid; }

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const { return _sizes.size(); }

       /*! \brief The number of elements in the mask. */
      size_type capacity() const { return _strides[dimension()]; }

      /*! \brief The lowest position in the grid. */
      const IndexArray& lower() const { return _lower; }

      /*! \brief The highest position in the grid. */
      const IndexArray& upper() const { return _upper; }

      /*! \brief Returns true if the set is empty. */
      bool empty() const { return std::find(_mask.begin(),_mask.end(),true)==_mask.end(); }

      /*! \brief The number of cells in the grid. */
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }

     /*! \brief A constant iterator to the beginning of the set. */
      const_iterator begin() const { return const_iterator(*this,0); }
      /*! \brief A constant iterator to the end of the set. */
      const_iterator end() const { return const_iterator(*this,capacity()); }

      /*! \brief Adjoins a cell to the set. */
      void adjoin(const GridCell<R>& c) {
        assert(c.grid()==this->_grid);
        compute_cell_mask(&_mask,_strides,_lower,c._position);
      }

      /*! \brief Adjoins a rectangle to the set. */
      void adjoin(const GridRectangle<R>& r) {
        assert(r.grid()==this->_grid);
        compute_rectangle_mask(&_mask,_strides,_lower,r._lower,r._upper);
      }

      /*! \brief Adjoins a GridMaskSet to the set. */
      void adjoin(const GridMaskSet<R>& ms) {
        assert(ms.grid()==this->_grid);
        _mask |= ms._mask;
      }

      /*! \brief Adjoins a GridCellListSet to the set. */
      void adjoin(const GridCellListSet<R>& cls) {
        assert(cls.grid()==this->_grid);
        compute_cell_list_mask(&_mask,_strides,_lower,cls._list);
      }

      /*! \brief Adjoins a GridRectangleListSet to the set. */
      void adjoin(const GridRectangleListSet<R>& rls) {
        assert(rls.grid()==this->_grid);
        compute_rectangle_list_mask(&_mask,_strides,_lower,rls._list);
      }

      friend std::ostream& operator<< <> (std::ostream&, const GridMaskSet<R>&);
      friend bool subset<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> join<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> regular_intersection<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
     private:
      const Grid<R>& _grid;
      array<index_type> _lower;
      array<index_type> _upper;
      array<size_type> _sizes;
      array<size_type> _strides;
      std::vector<bool> _mask;
    };



    /*! \brief A denotable set on a grid, defined using a list of cells.
     */
    template<typename R>
    class GridCellListSet {
      friend class GridRectangleListSet<R>;
      friend class GridMaskSet<R>;
     public:
      typedef R real_type;
      typedef State<R> state_type;

      /*!\brief Construct an empty set based on a Grid. */
      GridCellListSet(const Grid<R>& g);

      /*!\brief Construct from a GridMaskSet. */
      GridCellListSet(const GridMaskSet<R>& gms);

      /*!\brief Copy constructor. */
      GridCellListSet(const GridCellListSet<R>& g);

      /*!\brief Construct from a GridRectangleListSet. */
      GridCellListSet(const GridRectangleListSet<R>& g);

      /*!\brief Convert to a ListSet of Rectangles. */
      operator ListSet<R,Rectangle>() const;

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return _grid; }

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const { return _grid.dimension(); }

      /*! \brief True if the set is empty. */
      bool empty() const { return _list.empty(); }

      /*! \brief The numeber of cells in the list. */
      size_type size() const { return _list.size(); }

      /*!\brief The @a i th cell in the list. */
      GridCell<R> operator[] (const size_type i) const { return GridCell<R>(_grid,_list[i]); }

      /*!\brief Append a GridCell to the list. */
      void insert(const GridCell<R>& c) { _list.push_back(c._position); }

      friend std::ostream& operator<< <> (std::ostream&, const GridCellListSet<R>&);
     private:
      const Grid<R>& _grid;
      IntegerCellList _list;
    };





    /*! \brief A denotable set on a grid, defined using a list of rectangles.
     */
    template<typename R>
    class GridRectangleListSet {
      friend class GridMaskSet<R>;
      friend class GridCellListSet<R>;
     public:
      typedef R real_type;
      typedef State<R> state_type;

      /* TODO: Make iterators */

      /*!\brief Construct from a Grid. */
      GridRectangleListSet(const Grid<R>& g);

      /*!\brief Construct a set on a finer grid. */
      GridRectangleListSet(const Grid<R>& g, const GridCellListSet<R>& s);

      /*!\brief Construct a set on a finer grid. */
      GridRectangleListSet(const Grid<R>& g, const GridRectangleListSet<R>& s);

      /*!\brief Construct from a ListSet of Rectangles. */
      GridRectangleListSet(const ListSet<R,Rectangle>& s);

      /*!\brief Convert to a ListSet of Rectangles. */
      operator ListSet<R,Rectangle>() const;

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return _grid; }

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const { return _grid.dimension(); }

      /*! \brief True if the set is empty. */
      bool empty() const { return _list.empty(); }

      /*! \brief The number of rectangles in the list. */
      size_type size() const { return _list.size()/2; }

      /*!\brief Return the @a i th rectangle in the list. */
      GridRectangle<R> operator[] (const size_t i) const {
        return GridRectangle<R>(_grid,_list[2*i],_list[2*i+1]);
      }

      /*!\brief Append a GridRectangle to the list. */
      void insert(const GridRectangle<R>& r) {
        _list.push_back(r._lower);
        _list.push_back(r._upper);
      }

      friend std::ostream& operator<< <> (std::ostream&, const GridRectangleListSet<R>&);
     private:
      const Grid<R>& _grid;
      IntegerCellList _list;
    };



    /*! \brief Tests disjointness */
    template<typename R>
    bool disjoint (const GridMaskSet<R> &A, const GridMaskSet<R> &B);

    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect (const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return regular_intersection(A,B).empty();
    }

    /*! \brief Tests if A is a subset of the interior of B. */
    template<typename R>
    bool subset_of_interior (const GridMaskSet<R> &A, const GridMaskSet<R> &B);

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset (const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return (A._mask-B._mask).empty();
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    join(const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return GridMaskSet<R>(A._grid, A._lower, A._upper, A._mask | B._mask);
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    regular_intersection(const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return GridMaskSet<R>(A._grid, A._lower, A._upper, A._mask & B._mask);
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R> closure_of_intersection_of_interiors (const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      return regular_intersection(A,B);
    }




   /*! \brief A denotable set on a partition grid, defined using a partition tree of cells.
     */
    template<typename R>
    class PartitionTreeSet {
     public:
      /*! \brief Construct an empty set based on a PartitionGrid. */
      PartitionTreeSet(const PartitionGrid<R>& g)
        : _grid(g), _tree(1), _mask(1)
      {
        _tree[0]=PartitionGrid<R>::leaf;
        _mask[0]=false;
      }

      /*! \brief Construct a set based on a PartitionGrid, a subdivision tree and a mask. */
      PartitionTreeSet(const PartitionGrid<R>& g, const std::vector<bool>& t, const std::vector<bool>& m)
        : _grid(g), _tree(t), _mask(m)
      {
        assert(_tree.size() = 2*_mask.size()-1);
      }

      /*! \brief The space dimension of the set. */
      size_t dimension() const {
        return _grid.dimension();
      }
     private:
      const PartitionGrid<R>& _grid;
      std::vector<bool> _tree;
      std::vector<bool> _mask;
    };



    template<typename R>
    GridMaskSet_const_iterator<R>::GridMaskSet_const_iterator(const GridMaskSet<R>& s,
                                                              const IndexArray& pos)
      : _set(s),
        _index(compute_index(pos,s._lower,s._strides))
    {
      if(!_set._mask(_index)) {
        ++(*this);
      }
    }

    template<typename R>
    GridCell<R>
    GridMaskSet_const_iterator<R>::operator*() const
    {
      return GridCell<R>(_set.grid(),compute_position(_index,_set._lower,_set._strides));
    }

    template<typename R>
    GridMaskSet_const_iterator<R>&
    GridMaskSet_const_iterator<R>::operator++()
    {
      do {
        ++_index;
      } while( _index<_set._mask.size() && _set._mask[_index]==false );
      return *this;
    }





    template<typename R>
    GridMaskSet<R>::operator ListSet<R,Rectangle>() const
    {
      ListSet<R,Rectangle> result(dimension());
      for(typename GridMaskSet::const_iterator riter=begin(); riter!=end(); ++riter) {
        result.push_back(*riter);
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
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, IndexArray l, IndexArray u)
      : _grid(g), _lower(l), _upper(u),
        _sizes(g.dimension()), _strides(_sizes.size()+1),
        _mask()
    {
      _strides[0] = 1;
      for(dimension_type i=0; i!=g.dimension(); ++i) {
        assert(l[i] < u[i]);
        _sizes[i] = size_type(u[i] - l[i]);
        _strides[i+1] = _sizes[i] * _strides[i];
      }
      _mask.resize(_strides[g.dimension()]);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, IndexArray l, IndexArray u, const BooleanArray& m)
      : _grid(g), _lower(l), _upper(u),
        _sizes(g.dimension()), _strides(_sizes.size()+1),
        _mask(m)
    {
      _strides[0] = 1;
      for(dimension_type i=0; i!=g.dimension(); ++i) {
        assert(l[i] < u[i]);
        _sizes[i] = size_type(u[i] - l[i]);
        _strides[i+1] = _sizes[i] * _strides[i];
      }
      assert(m.size() == _strides[g.dimension()]);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
      : _grid(gcls.grid()), _lower(_grid.dimension()), _upper(_grid.dimension()),
        _sizes(_grid.dimension()), _strides(_grid.dimension()+1),
        _mask()
    {
      compute_cell_list_bounds(&_lower, &_upper, gcls._list);
      _sizes=_upper-_lower;
      _strides=compute_strides(_sizes);
      _mask=BooleanArray(_strides[dimension()],false);
      compute_cell_list_mask(&_mask, _strides, _lower, gcls._list);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridRectangleListSet<R>& grls)
      : _grid(grls.grid()), _lower(_grid.dimension()), _upper(_grid.dimension()),
        _sizes(_grid.dimension()), _strides(_grid.dimension()+1),
        _mask()
    {
      compute_rectangle_list_bounds(&_lower,&_upper, grls._list);
      _sizes=_upper-_lower;
      _strides=compute_strides(_sizes);
      _mask=BooleanArray(_strides[dimension()],false);
      compute_rectangle_list_mask(&_mask, _strides, _lower, grls._list);
    }



    template<typename R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
      : _grid(g), _list(g.dimension())
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
      : _grid(gms.grid()), _list(gms.dimension())
    {
      append_to_cell_list(&_list, gms._lower, gms._strides, gms._mask);
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridRectangleListSet<R>& grls)
      : _grid(grls.grid()), _list(grls.dimension())
    {
      append_to_cell_list(&_list, grls._list);
    }




    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g)
      : _grid(g), _list(g.dimension())
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridCellListSet<R>& s)
      : _grid(g), _list(g.dimension())
    {
      assert(g.dimension()==s.dimension());
      _list.resize(2*s.size());

      const FiniteGrid<R>* to=dynamic_cast<const FiniteGrid<R>*>(&g);
      const FiniteGrid<R>* from=dynamic_cast<const FiniteGrid<R>*>(&s._grid);

      array< std::vector<index_type> > tr = FiniteGrid<R>::index_translation(*from,*to);
      translate_cell_coordinates(&_list, s._list, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridRectangleListSet<R>& s)
      : _grid(g), _list(s._list)
    {
      const FiniteGrid<R>* to=dynamic_cast<const FiniteGrid<R>*>(&g);
      const FiniteGrid<R>* from=dynamic_cast<const FiniteGrid<R>*>(&s._grid);

      array< std::vector<index_type> > tr = FiniteGrid<R>::index_translation(*from,*to);
      translate_rectangle_coordinates(&_list, s._list, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const ListSet<R,Rectangle>& ls)
    /* FIXME: Memory leak! */
      : _grid(*new FiniteGrid<R>(ls)), _list(ls.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator ListSet_const_iterator;

      for(ListSet_const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        _list.push_back(_grid.lower_index(*riter));
        _list.push_back(_grid.upper_index(*riter));
      }
    }





    template <typename R>
    std::ostream& operator<<(std::ostream &os,
                             const GridRectangleListSet<R>& set)
    {
      os << "{ class: GridRectangleListSet<" << name<R>() << ">,\n  rectangle_list: [\n    ";
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
      os << "}" << std::endl;
      return os;
    }



    template <typename R>
    std::ostream& operator<<(std::ostream &os,
                             const GridCellListSet<R>& set)
    {
      os << "{ class: GridCellListSet<" << name<R>() << ">,\n  rectangle_list: [\n    ";
      if (!set.empty() ) {
        os << set[0];
      }
      for (size_t i=1; i<set.size(); ++i) {
        os << ",\n    " << Rectangle<R>(set[i]);
      }
      os << "\n  ]\n";
      os << "  grid: " << set.grid() << "\n";
      os << "  integer_cell_list: [ ";
      os.flush();
      for(size_type i=0; i!=set.size(); ++i) {
        if(i!=0) { os << ", "; }
        os <<  set[i].position();
      }
      os << " ]\n";
      os << "}" << std::endl;
      return os;
    }

    template <typename R>
    std::ostream& operator<<(std::ostream &os,
                             const GridMaskSet<R>& set)
    {
      os << "{ class: GridMaskSet<" << name<R>() << ">,\n";
      os << "  grid: " << set._grid << "\n";
      os << "  mask: " << set._mask << "\n";
      os << "}\n";
      return os;
    }

  }
}

#endif /* _GRID_DENOTABLE_SET_H */
