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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file grid_set.h
 *  \brief Denotable sets on grids.
 */

#ifndef _ARIADNE_GRID_SET_H
#define _ARIADNE_GRID_SET_H

#include <iosfwd>
#include <iostream>
#include <iterator>

#include "base/array.h"
#include "base/interval.h"
#include "base/binary_word.h"

#include "geometry/geometry_declarations.h"
#include "geometry/rectangle.h"
#include "geometry/point.h"
#include "geometry/grid_operations.h"

#include "geometry/list_set.h"

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

/*TODO: Unify bounds in FiniteGrid, and make GridMaskSet use only finite grid bounds*/

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

    template<typename R> bool interiors_intersect(const Rectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool interiors_intersect(const GridRectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const Rectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridRectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> regular_intersection(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> join(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> difference(const GridMaskSet<R>&, const GridMaskSet<R>&);

    template<typename R>
    GridRectangle<R>
    over_approximation(const Rectangle<R>& p, const Grid<R>& g);

    template<typename R>
    GridCellListSet<R>
    over_approximation(const Parallelopiped<R>& p, const Grid<R>& g);

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g); 

    template<typename R>
    GridRectangle<R>
    over_approximation_of_intersection(const Rectangle<R>& r1, 
                                       const Rectangle<R>& r2,
                                       const Grid<R>& g);
    
    template<typename R>
    GridCellListSet<R>
    over_approximation_of_intersection(const Parallelopiped<R>& p, 
                                       const Rectangle<R>& r,
                                       const Grid<R>& g);

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation_of_intersection(const ListSet<R,BS>& ls, 
                                       const Rectangle<R>& r, 
                                       const FiniteGrid<R>& g);
    
    template<typename R> std::ostream& operator<<(std::ostream&, const Grid<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridCell<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridRectangle<R>&);
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
      /*! \brief The interval @math [p_n,p_{n+1}) in dimension @a d index containing @a x. */
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
          res[i]=subdivision_lower_index(i,r.lower_bound(i));
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


    /*!\brief A finite, nonuniform grid of rectangles in Euclidean space. */
    template<typename R>
    class FiniteGrid : public Grid<R> {
      typedef R real_type;
     public:
      /*! \brief Cloning operator. */
      virtual FiniteGrid<R>* clone() const;
    
      /*! \brief Destructor. */
      virtual ~FiniteGrid();

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { 
        return this->_subdivision_coordinates.size();
      }
    
      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const {
        return _subdivision_coordinates[d][n];
      }
  
      /*! \brief The index of interval in dimension @a d index containing @a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        typename std::vector<R>::const_iterator pos;
        if(x<_subdivision_coordinates[d].front() || x>_subdivision_coordinates[d].back()) {
          throw std::runtime_error("point does not lie in extent of finite grid");
        }
        pos = std::upper_bound(_subdivision_coordinates[d].begin(),
                               _subdivision_coordinates[d].end(), x);
        return (pos - _subdivision_coordinates[d].begin()) - 1;
      }


      /*! \brief Construct from a bounding box and an equal number of subdivisions in each coordinate. */
      explicit FiniteGrid(const Rectangle<R>& r, size_type n);

      /*! \brief Construct from a bounding box and an array giving the number of subdivisions in each coordinate. */
      explicit FiniteGrid(const Rectangle<R>& r, SizeArray sz);

      /*! \brief Construct from a list of subdivision coordinates in each dimension. */
      explicit FiniteGrid(const array< std::vector<R> >& sp);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit FiniteGrid(const ListSet<R,Rectangle>& ls);

      /*! \brief Join two finite grids. */
      FiniteGrid(const FiniteGrid& g1, FiniteGrid& g2);

      bool operator==(const FiniteGrid<R>& g) const { 
        return this->_subdivision_coordinates==g._subdivision_coordinates; }
            
      /*! \brief The total number of cells. */
      size_type capacity() const { return _strides[dimension()]; }
      /*! \brief The number of subdivision intervals in dimension @a d. */
      size_type size(dimension_type d) const { return _subdivision_coordinates[d].size()-1; }

      /*! \brief The lowest valid vertex index. */
      IndexArray lower() const { return IndexArray(dimension(),0); }
      /*! \brief The highers valid vertex index. */
      IndexArray upper() const { return lower()+sizes(); }
      /*! \brief The block of valid lattice cells. */
      IndexBlock bounds() const { return IndexBlock(lower(),upper()); }
      /*! \brief The number of subdivision intervals in each dimension. */
      SizeArray sizes() const {
        SizeArray result(dimension());
        for(dimension_type d=0; d!=dimension(); ++d) {
          result[d]=_subdivision_coordinates[d].size()-1;
        }
        return result;
      }

      Rectangle<R> bounding_box() { 
        Rectangle<R> result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result.set_lower_bound(i,this->_subdivision_coordinates[i].front());
          result.set_upper_bound(i,this->_subdivision_coordinates[i].back());
        }
        return result;
      }
      
      /*! \brief Find the rule to translate elements from a grid to a refinement. */
      static array< std::vector<index_type> > index_translation(const FiniteGrid<R>& from, const FiniteGrid<R>& to);

      virtual std::ostream& write(std::ostream& os) const;
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
      InfiniteGrid<R>* clone() const;
      
      /*! \brief Construct from an array of subdivision lengths \a sl.
       */
      InfiniteGrid(const array<R>& sl)
        : _subdivision_lengths(sl) { }

      /*! \brief The underlying dimension of the grid. */
      virtual dimension_type dimension() const { 
        return _subdivision_lengths.size(); }

      /*! \brief The coordinate of the @a n th subdivision point in dimension @a d. */
      virtual real_type subdivision_coordinate(dimension_type d, index_type n) const { 
        return _subdivision_lengths[d] * n; }

      /*! \brief The length of the subdivision in the \a d th coordinate. */
      real_type subdivision_length(dimension_type d) const { 
        return _subdivision_lengths[d]; }

      /*! \brief The index of interval in dimension @a d index containing @a x. */
      virtual index_type subdivision_interval(dimension_type d, const real_type& x) const {
        index_type result = quotient(x,_subdivision_lengths[d]);
        return result;
      }

      virtual std::ostream& write(std::ostream& os) const;
     private:
      array<R> _subdivision_lengths;
    };



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

      /*!\brief Construct an empty rectangle on a grid. */
      GridRectangle(const Grid<R>& g);
      /*!\brief Construct from a grid and a bounding block. */
      GridRectangle(const Grid<R>& g, const IndexBlock& b);
      /*!\brief Construct from a grid and two integer arrays giving the corners. */
      GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u);
      /*!\brief Construct from a grid and an ordinary rectangle. */
      GridRectangle(const Grid<R>& g, const Rectangle<R>& r);
      /*!\brief Construct from a GridCell. */
      GridRectangle(const GridCell<R>& c);

      /*!\brief The grid containing the rectangle. */
      const Grid<R>& grid() const { return _grid; }
      /*!\brief The dimension of the rectangle. */
      dimension_type dimension() const { return _position.dimension(); }

      /*!\brief The position of the rectangle in the grid. */
      const IndexBlock& position() const { return _position; }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Grid<R>& _grid;
      IndexBlock _position;
    };


    template<class R>
    class GridMaskSetConstIterator {
      typedef GridMaskSetConstIterator Self;
  
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridCell<R> value_type;
      typedef GridCell<R> reference;
      typedef const GridCell<R>* pointer;
      typedef int difference_type;
     public:
      GridMaskSetConstIterator(const GridMaskSet<R>& s, const IndexArray& pos);
      GridMaskSetConstIterator(const GridMaskSet<R>& s, const size_type& ind);
      GridCell<R> operator*() const;
      Self& operator++();
      Self operator++(int) { Self tmp=*this; ++(*this); return tmp; }
      bool operator==(const Self& other) const { return _index==other._index && (&_set)==(&(other._set)); }
      bool operator!=(const Self& other) const { return !(*this==other); }
    private:
      void initialise();
    public:
      const GridMaskSet<R>& _set;
      size_type _index;
    };

    
    
    /*! \brief A denotable set on a finite grid, defined using a mask. */
    template<typename R>
    class GridMaskSet {
      friend class GridMaskSetConstIterator<R>;
      friend class GridCellListSet<R>;
     public:
      typedef R real_type;
      typedef Point<R> state_type;
      typedef GridMaskSetConstIterator<R> iterator;
      typedef GridMaskSetConstIterator<R> const_iterator;

      /*!\brief Construct an empty set from a finite grid. */
      GridMaskSet(const FiniteGrid<R>& g);
     
      /*!\brief Construct an empty set from a grid, and a bounding box. */
      GridMaskSet(const Grid<R>& g, const IndexBlock& b);
     
      /*!\brief Construct a set from a grid, a bounding box, and a mask */
      GridMaskSet(const Grid<R>& g, const IndexBlock& b, const BooleanArray& m);

      /*!\brief Copy constructor. */
      GridMaskSet(const GridMaskSet<R>& gms);

      /*!\brief Construct from a ListSet of Rectangle. */
      GridMaskSet(const ListSet<R,Rectangle>& ls);

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
        if(this->grid()==gms.grid() && this->_bounds==gms._bounds) {
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
      const IndexArray& lower() const { return _bounds.lower(); }

      /*! \brief The highest position in the grid. */
      const IndexArray& upper() const { return _bounds.upper(); }

      /*! \brief The highest position in the grid. */
      const IndexBlock& bounds() const { return _bounds; }

      /*! \brief The number of cells in the grid. */
      const BooleanArray& mask() const { return _mask; }

      /*! \brief Returns true if the set is empty. */
      bool empty() const { return std::find(_mask.begin(),_mask.end(),true)==_mask.end(); }

      /*! \brief The number of cells in the grid. */
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }

      /*! \brief The ith nonempty cell in the grid. */
      GridCell<R> operator[](size_type i) const { 
        const_iterator iter=this->begin();
        while(i!=0) {
          --i;
          ++iter;
        }
        return *iter;
      }

      /*! \brief A constant iterator to the beginning of the set. */
      const_iterator begin() const { return const_iterator(*this,0); }
      /*! \brief A constant iterator to the end of the set. */
      const_iterator end() const { return const_iterator(*this,capacity()); }

      Rectangle<R> bounding_box() const {
        return GridRectangle<R>(grid(),bounds());
      }
      
      /*! \brief Adjoins a cell to the set. */
      void adjoin(const GridCell<R>& c) {
        assert(c.grid()==this->grid());
        compute_cell_mask(&_mask,_strides,_lower,c._position);
      }

      /*! \brief Adjoins a rectangle to the set. */
      void adjoin(const GridRectangle<R>& r) {
        assert(r.grid()==this->grid());
        compute_rectangle_mask(&_mask,_strides,_lower,
                               r.position().lower(),r.position().upper());
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

      friend bool subset<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> join<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> regular_intersection<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> difference<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      IndexBlock _bounds;
      IndexArray _lower;
      SizeArray _sizes;
      SizeArray _strides;
      BooleanArray _mask;
      //size_type _capacity;
      //index_type _offset;
    };



    template<class R>
    class GridCellListSetConstIterator {
      typedef GridCellListSetConstIterator Self;
     public:
      typedef std::random_access_iterator_tag iterator_category;
      typedef GridCell<R> value_type;
      typedef GridCell<R> reference;
      typedef const GridCell<R>* pointer;
      typedef int difference_type;
     public:
      GridCellListSetConstIterator(const GridCellListSet<R>& s, const size_type& ind) : _set(s), _index(ind) { }
      GridCell<R> operator*() const { return _set[_index]; }
      Self& operator++() { ++_index; return *this; }
      Self operator++(int) { Self tmp=*this; ++(*this); return tmp; }
      bool operator==(const Self& other) const { return _index==other._index && (&_set)==(&(other._set)); }
      bool operator!=(const Self& other) const { return !(*this==other); }
    public:
      const GridCellListSet<R>& _set;
      size_type _index;
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
      typedef GridCellListSetConstIterator<R> iterator;
      typedef GridCellListSetConstIterator<R> const_iterator;

      /*!\brief Construct an empty set based on a Grid. */
      GridCellListSet(const Grid<R>& g);

      /*!\brief Construct from a GridMaskSet. */
      GridCellListSet(const GridMaskSet<R>& gms);

      /*!\brief Copy constructor. */
      GridCellListSet(const GridCellListSet<R>& gcls);

      /*!\brief Construct from a GridRectangleListSet. */
      GridCellListSet(const GridRectangleListSet<R>& grls);

      /*!\brief Construct from a GridRectangleListSet. */
      GridCellListSet(const ListSet<R,Rectangle>& rls);

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

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const { return const_iterator(*this,0); }

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const { return const_iterator(*this,this->size()); }

      /*!\brief Append a GridCell to the list. */
      void adjoin(const GridCell<R>& c) { _list.push_back(c._position); }

      friend std::ostream& operator<< <> (std::ostream&, const GridCellListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      IndexArrayList _list;
    };





    template<class R>
    class GridRectangleListSetConstIterator {
      typedef GridRectangleListSetConstIterator Self;
     public:
      typedef std::random_access_iterator_tag iterator_category;
      typedef GridRectangle<R> value_type;
      typedef GridRectangle<R> reference;
      typedef const GridRectangle<R>* pointer;
      typedef int difference_type;
     public:
      GridRectangleListSetConstIterator(const GridRectangleListSet<R>& s, const size_type& ind) : _set(s), _index(ind) { }
      GridRectangle<R> operator*() const { return _set[_index]; }
      Self& operator++() { ++_index; return *this; }
      Self operator++(int) { Self tmp=*this; ++(*this); return tmp; }
      bool operator==(const Self& other) const { return _index==other._index && (&_set)==(&(other._set)); }
      bool operator!=(const Self& other) const { return !(*this==other); }
    public:
      const GridRectangleListSet<R>& _set;
      size_type _index;
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
      typedef GridRectangleListSetConstIterator<R> iterator;
      typedef GridRectangleListSetConstIterator<R> const_iterator;

      /*!\brief Destructor. */
      ~GridRectangleListSet() { 
        //delete _grid_ptr; }
      }
      
      /*!\brief Construct from a Grid. */
      GridRectangleListSet(const Grid<R>& g);

      /*!\brief Construct from a ListSet of Rectangles. */
      explicit GridRectangleListSet(const ListSet<R,Rectangle>& s);

      /*!\brief Construct from a PartitionTreeSet. */
      /* FIXME: This constructor is only included since Boost Python doesn't find the conversion operator */
      explicit GridRectangleListSet(const PartitionTreeSet<R>& s);

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
        return GridRectangle<R>(grid(),IndexBlock(_list[2*i],_list[2*i+1]));
      }

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const { return const_iterator(*this,0); }

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const { return const_iterator(*this,this->size()); }

      /*!\brief Append a GridRectangle to the list. */
      void push_back(const GridRectangle<R>& r) {
        _list.push_back(r._position.lower());
        _list.push_back(r._position.upper());
      }

      friend std::ostream& operator<< <> (std::ostream&, const GridRectangleListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      IndexBlockList _list;
    };


  }
}

#endif /* _GRID_SET_H */
