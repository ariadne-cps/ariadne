/***************************************************************************
 *            grid_set.h
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
 
/*! \file grid_set.h
 *  \brief Denotable sets on grids.
 */

#ifndef _ARIADNE_GRID_SET_H
#define _ARIADNE_GRID_SET_H

#include <iostream>

#include "array.h"
#include "interval.h"
#include "binary_word.h"
#include "rectangle.h"
#include "point.h"

#include "grid_operations.h"

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
    typedef Rectangle<index_type> IntegerRectangle;
    
    template<typename R, template<typename> class BS> class ListSet;

    template<typename R> class Rectangle;
    template<typename R> class Point;

    template<typename R> class Grid;
    template<typename R> class FiniteGrid;
    template<typename R> class UniformGrid;

    template<typename R> class GridCell;
    template<typename R> class GridRectangle;
    template<typename R> class GridMaskSet;
    template<typename R> class GridRectangleListSet;
    template<typename R> class GridCellListSet;

    template<typename R> bool subset(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> regular_intersection(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> join(const GridMaskSet<R>&, const GridMaskSet<R>&);

    template<typename R> std::ostream& operator<<(std::ostream&, const GridRectangleListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridCellListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridMaskSet<R>&);
    
    /*!\brief Base type for defining a grid.
     * We use inheritence and abstract functions here partly for ease of development
     * and partly since the grid coordinates only play a role when converting to rectangles
     * and should occur with linear complexity in the space dimension.
     */
    template<typename R>
    class Grid {
    public:
      typedef R real_type;
     
      /*! \brief Destructor. */
      virtual ~Grid() { }

      /*! \brief Cloning operator. */
      virtual Grid<R>* clone() const = 0;
      
      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const = 0;
      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const = 0;
      /*! \brief The interval @math [p_n,p_{n+1}] in dimension @a d index containing @a x. */
     virtual index_type subdivision_interval(dimension_type d, const real_type& x) const = 0;

      /*! \brief The index in dimension @a d of the subdivision point @a x. Throws an exception if @a x is not a subdivision point. */
      index_type subdivision_index(dimension_type d, const real_type& x) const {
        index_type n=subdivision_interval(d,x);
        if(subdivision_coordinate(d,n) == x) { return n; }
        throw std::domain_error("Value is not a grid coordinate");
      }

      /*! \brief The index of the subdivision point below x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const {
        return subdivision_interval(d,x);
      }

      /*! \brief The index of the subdivision point above x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const {
        index_type n=subdivision_interval(d,x);
        return subdivision_coordinate(d,n) == x ? n : n+1;
      }

      bool operator==(const Grid<R>& g) const { return this==&g; }

      IndexArray index(const Point<R>& s) const {
        IndexArray res(s.dimension());
        for(size_type i=0; i!=res.size(); ++i) {
          res[i]=subdivision_index(i,s[i]);
        }
        return res;
      }


      IndexArray lower_index(const Rectangle<R>& r) const {
        IndexArray res(r.dimension());
        for(size_type i=0; i!=res.size(); ++i) {
          res[i]=subdivision_lower_index(i,r.lower(i));
        }
        return res;
      }

      IndexArray upper_index(const Rectangle<R>& r) const {
        IndexArray res(r.dimension());
        for(size_type i=0; i!=res.size(); ++i) {
          res[i]=subdivision_upper_index(i,r.upper_bound(i));
        }
        return res;
      }

      virtual std::ostream& write(std::ostream& os) const = 0;
   };


    /*!\brief A finite, nonuniform grid of rectangles in Euclidean space.
     */
    template<typename R>
    class FiniteGrid : public Grid<R> {
      typedef R real_type;
     public:
      /*! \brief Cloning operator. */
      FiniteGrid<R>* clone() const { return new FiniteGrid<R>(*this); }
      
      /*! \brief Destructor. */
      virtual ~FiniteGrid() { }

      /*! \brief Construct from a list of subdivision coordinates in each dimension. */
      explicit FiniteGrid(const array< std::vector<R> >& sp);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit FiniteGrid(const ListSet<R,Rectangle>& ls);

      /*! \brief Join two finite grids. */
      FiniteGrid(const FiniteGrid& g1, FiniteGrid& g2);

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { return _subdivision_coordinates.size(); }
      /*! \brief The total number of cells. */
      size_type capacity() const { return _strides[dimension()]; }
      /*! \brief The number of subdivision intervals in dimension @a d. */
      size_type size(dimension_type d) const { return _subdivision_coordinates[d].size(); }

      /*! \brief The lowest valid vertex index. */
      IndexArray lower() const { return IndexArray(dimension(),0); }
      /*! \brief The highers valid vertex index. */
      IndexArray upper() const { return lower()+sizes(); }
      /*! \brief The number of subdivision intervals in each dimension. */
      SizeArray sizes() const {
        SizeArray result(dimension());
        for(dimension_type d=0; d!=dimension(); ++d) {
          result[d]=_subdivision_coordinates[d].size();
        }
        return result;
      }

      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const {
        return _subdivision_coordinates[d][n];
      }
      /*! \brief The index of interval in dimension @a d index containing @a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        typename std::vector<R>::const_iterator pos;
        pos = std::upper_bound(_subdivision_coordinates[d].begin(),
                               _subdivision_coordinates[d].end(), x);
        return (pos - _subdivision_coordinates[d].begin()) - 1;
      }

      /*! \brief Find the rule to translate elements from a grid to a refinement. */
      static array< std::vector<index_type> > index_translation(const FiniteGrid<R>& from, const FiniteGrid<R>& to);

      virtual std::ostream& write(std::ostream& os) const {
        return os << *this;
      }
     private:
      void create();
     private:
      array< std::vector<R> > _subdivision_coordinates;
      SizeArray _strides;
    };



    /*!\brief An infinite, uniform grid of rectangles in Euclidean space.
     */
    template<typename R> class InfiniteGrid : public Grid<R> {
      typedef R real_type;
     public:
      /*! \brief Cloning operator. */
      InfiniteGrid<R>* clone() const { return new InfiniteGrid<R>(*this); }
      
      /*! \brief Construct from an array of subdivision lengths \a sl.
       */
      InfiniteGrid(const array<R>& sl)
        : _subdivision_lengths(sl) { }

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { return _subdivision_lengths.size(); }

//      const array<R>& subdivision_lengths() const { return _subdivision_lengths; }
      /*! \brief The length of the subdivision in the \a d th coordinate. */
      real_type subdivision_length(dimension_type d) const { return _subdivision_lengths[d]; }
      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const { return _subdivision_lengths[d] * n; }
      virtual real_type subdivision_coordinate(dimension_type d, size_type n) const { return _subdivision_lengths[d] * n; }
      /*! \brief The index of interval in dimension @a d index containing @a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        index_type result = quotient(x,_subdivision_lengths[d]);
        return result;
      }
      virtual std::ostream& write(std::ostream& os) const {
        return os << *this;
      }
     private:
      array<R> _subdivision_lengths;
    };

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
    operator<<(std::ostream& os, const FiniteGrid<R>& g)
    {
      os << "[ ";
      for(dimension_type n=0; n!=g.dimension(); ++n) {
        if(n!=0) { os << " x "; }
        os << "[";
        for(size_type i=0; i!=g.size(n); ++i) {
          if(i!=0) { os << ", "; }
          os << g.subdivision_coordinate(n,i);
        }
        os << "]";
      }
      os << " ]";
      return os;
    }

    /*! \brief A cell in a grid.
     */
    template<typename R>
    class GridCell {
      friend class GridRectangle<R>;
      friend class GridMaskSet<R>;
      friend class GridCellListSet<R>;
     public:
      typedef R real_type;
      typedef Point<R> state_type;

      /*!\brief Construct from a grid and an integer array. */
      GridCell(const Grid<R>& g, const IndexArray& pos);

      /*!\brief The grid containing the cell. */
      const Grid<R>& grid() const { return _grid; }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return _position.size(); }

      /*!\brief The position of the cell in the grid. */
      const IndexArray& position() const { return _position; }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Grid<R>& _grid;
      IndexArray _position;
    };


    /*! \brief A rectangle in a grid.
     */
    template<typename R>
    class GridRectangle {
      friend class GridCell<R>;
      friend class GridMaskSet<R>;
      friend class GridRectangleListSet<R>;
     public:
      typedef R real_type;
      typedef Point<R> state_type;

      /*!\brief Construct from a grid and two integer arrays giving the corners. */
      GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u);
      /*!\brief Construct from a grid and an ordinary rectangle. */
      GridRectangle(const Grid<R>& g, const Rectangle<R>& r);
      /*!\brief Construct from a GridCell. */
      GridRectangle(const GridCell<R>& c);

      /*!\brief The grid containing the rectangle. */
      const Grid<R>& grid() const { return _grid; }
      /*!\brief The dimension of the rectangle. */
      dimension_type dimension() const { return _lower.size(); }

      /*!\brief The position of the rectangle in the grid. */
      IntegerRectangle position() const {
        Point<index_type> lst(_lower.begin(),_lower.end());
        Point<index_type> ust(_upper.begin(),_upper.end());
        return IntegerRectangle(lst,ust);
      }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Grid<R>& _grid;
      array<index_type> _lower;
      array<index_type> _upper;
    };


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
      typedef Point<R> state_type;
      typedef GridMaskSet_const_iterator<R> const_iterator;

      /*!\brief Construct an empty set from a grid, and a list of lower and upper bounds. */
      GridMaskSet(const Grid<R>& g, IndexArray l, IndexArray u);

      /*!\brief Construct an empty set from a grid, a list of lower and upper bounds, and a mask */
      GridMaskSet(const Grid<R>& g, IndexArray l, IndexArray u, const BooleanArray& m);

      /*!\brief Copy constructor. */
      GridMaskSet(const GridMaskSet<R>& gms);

      /*!\brief Construct from a GridCellListSet. */
      GridMaskSet(const GridCellListSet<R>& gcls);

      /*!\brief Construct from a GridRectangleListSet. */
      GridMaskSet(const GridRectangleListSet<R>& grls);

      /*!\brief Construct from a Grid and a ListSet of Rectangle. */
      GridMaskSet(const FiniteGrid<R>&, const ListSet<R,Rectangle>& ls);

      /*!\brief Convert to a ListSet of Rectangles. */
      operator ListSet<R,Rectangle>() const;

      /*! \brief Equality operator. */
      bool operator==(const GridMaskSet<R>& gms) {
        if(this->grid()==gms.grid() && this->_lower==gms._lower && this->_upper==gms._upper) {
          return this->_mask==gms._mask;
        }
        throw(std::domain_error("Can only compare GridMaskSets on the same grid"));
      }

      /*! \brief Inequality operator. */
      bool operator!=(const GridMaskSet<R>& gms) { return !(*this==gms); }

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return *_grid_ptr; }

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

       /*! \brief Adjoins a rectangle to the set. */
      void adjoin(const Rectangle<R>& r) {
        IndexArray rlower=grid().index(r.lower_corner());
        IndexArray rupper=grid().index(r.upper_corner());
        compute_rectangle_mask(&_mask,_strides,_lower,rlower,rupper);
      }

      /*! \brief Adjoins a cell to the set. */
      void adjoin(const GridCell<R>& c) {
        assert(c.grid()==this->grid());
        compute_cell_mask(&_mask,_strides,_lower,c._position);
      }

      /*! \brief Adjoins a rectangle to the set. */
      void adjoin(const GridRectangle<R>& r) {
        assert(r.grid()==this->grid());
        compute_rectangle_mask(&_mask,_strides,_lower,r._lower,r._upper);
      }

      /*! \brief Adjoins a GridMaskSet to the set. */
      void adjoin(const GridMaskSet<R>& ms) {
        assert(ms.grid()==this->grid());
        _mask |= ms._mask;
      }

      /*! \brief Adjoins a GridCellListSet to the set. */
      void adjoin(const GridCellListSet<R>& cls) {
        assert(cls.grid()==this->grid());
        compute_cell_list_mask(&_mask,_strides,_lower,cls._list);
      }

      /*! \brief Adjoins a GridRectangleListSet to the set. */
      void adjoin(const GridRectangleListSet<R>& rls) {
        assert(rls.grid()==this->grid());
        compute_rectangle_list_mask(&_mask,_strides,_lower,rls._list);
      }

      friend std::ostream& operator<< <> (std::ostream&, const GridMaskSet<R>&);
      friend bool subset<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> join<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> regular_intersection<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
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
      typedef Point<R> state_type;

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
      const Grid<R>& grid() const { return *_grid_ptr; }

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const { return grid().dimension(); }

      /*! \brief True if the set is empty. */
      bool empty() const { return _list.empty(); }

      /*! \brief The numeber of cells in the list. */
      size_type size() const { return _list.size(); }

      /*!\brief The @a i th cell in the list. */
      GridCell<R> operator[] (const size_type i) const { return GridCell<R>(grid(),_list[i]); }

      /*!\brief Append a GridCell to the list. */
      void insert(const GridCell<R>& c) { _list.push_back(c._position); }

      friend std::ostream& operator<< <> (std::ostream&, const GridCellListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
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
      typedef Point<R> state_type;

      /* TODO: Make iterators */

      /*!\brief Destructor. */
      ~GridRectangleListSet() { delete _grid_ptr; }

      /*!\brief Construct from a Grid. */
      GridRectangleListSet(const Grid<R>& g);

      /*!\brief Construct from a ListSet of Rectangles. */
      GridRectangleListSet(const ListSet<R,Rectangle>& s);

      /*!\brief Copy constructor. */
      GridRectangleListSet(const GridRectangleListSet<R>& s);

      /*!\brief Construct a set on a finer grid. */
      GridRectangleListSet(const Grid<R>& g, const ListSet<R,Rectangle>& s);

      /*!\brief Construct a set on a finer grid. */
      GridRectangleListSet(const Grid<R>& g, const GridCellListSet<R>& s);

      /*!\brief Construct a set on a finer grid. */
      GridRectangleListSet(const Grid<R>& g, const GridRectangleListSet<R>& s);

      /*!\brief Convert to a ListSet of Rectangles. */
      operator ListSet<R,Rectangle>() const;

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return *_grid_ptr; }

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const { return grid().dimension(); }

      /*! \brief True if the set is empty. */
      bool empty() const { return _list.empty(); }

      /*! \brief The number of rectangles in the list. */
      size_type size() const { return _list.size()/2; }

      /*!\brief Return the @a i th rectangle in the list. */
      GridRectangle<R> operator[] (const size_t i) const {
        return GridRectangle<R>(grid(),_list[2*i],_list[2*i+1]);
      }

      /*!\brief Append a GridRectangle to the list. */
      void adjoin(const GridRectangle<R>& r) {
        _list.push_back(r._lower);
        _list.push_back(r._upper);
      }

      friend std::ostream& operator<< <> (std::ostream&, const GridRectangleListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      IntegerCellList _list;
    };



    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid(g), _lower(l), _upper(u)
    {
      assert(g.dimension()==l.size());
      assert(l.size()==u.size());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const Rectangle<R>& r)
      : _grid(g), _lower(g.dimension()), _upper(g.dimension())
    {
      assert(g.dimension()==r.dimension());
      for(dimension_type d=0; d!=dimension(); ++d) {
        /* TODO: Catch and rethrow exceptions */
        _lower[d]=g.subdivision_index(d,r.lower_bound(d));
        _upper[d]=g.subdivision_index(d,r.upper_bound(d));
      }
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const GridCell<R>& c)
      : _grid(c._grid), _lower(c._position), _upper(c._position)
    {
      for(dimension_type i=0; i!=_upper.size(); ++i) {
        _upper[i]+=1;
      }
    }


    template<typename R>
    GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
      : _grid(g), _position(pos)
    {
      assert(_position.size()==_grid.dimension());
    }



    template<typename R>
    inline
    GridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(dimension_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_position[i]));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_position[i]+1));
      }

      return result;
    }


    template<typename R>
    inline
    GridRectangle<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(size_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_lower[i]));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_upper[i]));
      }

      return result;
    }




    template<typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const GridRectangle<R>& r) {
      return os << Rectangle<R>(r);
    }

    template<typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const GridCell<R>& c) {
      return os << Rectangle<R>(c);
    }



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

    /*! \brief The union of \a A and \a B. */
    template<typename R>
    GridMaskSet<R>
    join(const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return GridMaskSet<R>(A.grid(), A._lower, A._upper, A._mask | B._mask);
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    regular_intersection(const GridMaskSet<R> &A, const GridMaskSet<R> &B)
    {
      assert(A.grid()==B.grid());
      return GridMaskSet<R>(A.grid(), A._lower, A._upper, A._mask & B._mask);
    }




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
      : _grid_ptr(g.clone()), _lower(l), _upper(u),
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
      : _grid_ptr(g.clone()), _lower(l), _upper(u),
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
    GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
      : _grid_ptr(gms.grid().clone()), _lower(gms._lower), _upper(gms._upper),
        _sizes(gms._sizes), _strides(gms._strides),
        _mask(gms._mask)
    {  
    }
    
    
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(gcls.grid().clone()), _lower(grid().dimension()), _upper(grid().dimension()),
        _sizes(grid().dimension()), _strides(grid().dimension()+1),
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
      : _grid_ptr(grls.grid().clone()), _lower(grid().dimension()), _upper(grid().dimension()),
        _sizes(grid().dimension()), _strides(grid().dimension()+1),
        _mask()
    {
      compute_rectangle_list_bounds(&_lower,&_upper, grls._list);
      _sizes=_upper-_lower;
      _strides=compute_strides(_sizes);
      _mask=BooleanArray(_strides[dimension()],false);
      compute_rectangle_list_mask(&_mask, _strides, _lower, grls._list);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg, const ListSet<R,Rectangle>& ls)
      : _grid_ptr(fg.clone()), _lower(fg.lower()), _upper(fg.upper()),
        _sizes(fg.dimension()), _strides(fg.dimension()+1),
        _mask(fg.capacity())
    {
      _sizes=_upper-_lower;
      _strides=compute_strides(_sizes);
      _mask=BooleanArray(_strides[dimension()],false);
      for(typename ListSet<R,Rectangle>::const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        GridRectangle<R> r(grid(),*riter);
        adjoin(r);
      }
    }



    template<typename R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
      : _grid_ptr(g.clone()), _list(g.dimension())
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
      : _grid_ptr(gms.grid().clone()), _list(gms.dimension())
    {
      append_to_cell_list(&_list, gms._lower, gms._strides, gms._mask);
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridRectangleListSet<R>& grls)
      : _grid_ptr(grls.grid().clone()), _list(grls.dimension())
    {
      append_to_cell_list(&_list, grls._list);
    }




    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g)
      : _grid_ptr(g.clone()), _list(g.dimension())
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const ListSet<R,Rectangle>& s)
      : _grid_ptr(new FiniteGrid<R>(s)), _list(s.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator list_set_const_iterator;
      
      _list.resize(2*s.size());

      for(list_set_const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
        this->adjoin(GridRectangle<R>(grid(),*iter));
      }
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const GridRectangleListSet<R>& s)
      : _grid_ptr(s.grid().clone()), _list(s._list)
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridCellListSet<R>& s)
      : _grid_ptr(g.clone()), _list(g.dimension())
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
      : _grid_ptr(g.clone()), _list(s._list)
    {
      const FiniteGrid<R>* to=dynamic_cast<const FiniteGrid<R>*>(&g);
      const FiniteGrid<R>* from=dynamic_cast<const FiniteGrid<R>*>(&s.grid());

      array< std::vector<index_type> > tr = FiniteGrid<R>::index_translation(*from,*to);
      translate_rectangle_coordinates(&_list, s._list, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const ListSet<R,Rectangle>& ls)
      : _grid_ptr(g.clone()), _list(ls.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator ListSet_const_iterator;

      for(ListSet_const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        _list.push_back(grid().lower_index(*riter));
        _list.push_back(grid().upper_index(*riter));
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
      os << "  grid: " << set.grid() << "\n";
      os << "  mask: " << set._mask << "\n";
      os << "}\n";
      return os;
    }

  }
}

#endif /* _GRID_SET_H */
