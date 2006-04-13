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
 *over
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

#include <boost/iterator/iterator_adaptor.hpp>

#include "../declarations.h"

#include "../base/array.h"

#include "../geometry/rectangle.h"
#include "../geometry/grid.h"
#include "../geometry/lattice_set.h"

/*TODO: Unify bounds in FiniteGrid, and make GridMaskSet use only finite grid bounds*/

namespace Ariadne {
  namespace Geometry {

    template<typename R> bool interiors_intersect(const GridRectangle<R>&, const GridRectangle<R>&);
    template<typename R> bool interiors_intersect(const Rectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool interiors_intersect(const GridRectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool interiors_intersect(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridRectangle<R>&, const GridRectangle<R>&);
    template<typename R> bool subset(const Rectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridRectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool inner_subset(const Rectangle<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> regular_intersection(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> join(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> difference(const GridMaskSet<R>&, const GridMaskSet<R>&);

    template<typename R>
    GridRectangle<R>
    outer_approximation(const Rectangle<R>& p, const Grid<R>& g);

    template<typename R>
    GridRectangle<R>
    over_approximation(const Rectangle<R>& p, const Grid<R>& g);

    template<typename R>
    GridCellListSet<R>
    over_approximation(const Parallelotope<R>& p, const Grid<R>& g);

    template<typename R>
    GridCellListSet<R>
    over_approximation(const Zonotope<R>& p, const Grid<R>& g);

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
    over_approximation_of_intersection(const Parallelotope<R>& p, 
                                       const Rectangle<R>& r,
                                       const Grid<R>& g);

    template<typename R>
    GridCellListSet<R>
    over_approximation_of_intersection(const Zonotope<R>& p, 
                                       const Rectangle<R>& r,
                                       const Grid<R>& g);

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation_of_intersection(const ListSet<R,BS>& ls, 
                                       const Rectangle<R>& r, 
                                       const FiniteGrid<R>& g);
    
    template<typename R> std::ostream& operator<<(std::ostream&, const GridCell<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridRectangle<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridRectangleListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridCellListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridMaskSet<R>&);
    
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
      
      /*!\brief Construct from a grid and an unit grid cell. */
      GridCell(const Grid<R>& g, const LatticeCell& pos);

      /*!\brief Construct from a grid and an unit grid cell. */
      GridCell(const Grid<R>& g, const IndexArray& pos);

      /*!\brief The grid containing the cell. */
      const Grid<R>& grid() const { return _grid; }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return _position.dimension(); }

      /*!\brief The position of the cell in the grid. */
      const LatticeCell& position() const { return _position; }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Grid<R>& _grid;
      LatticeCell _position;
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
      GridRectangle(const Grid<R>& g, const LatticeRectangle& b);
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
      const LatticeRectangle& position() const { return _position; }

      /*!\brief Tests if the rectangle is empty. */
      bool empty() const { return _position.empty(); }
      /*!\brief Tests if the rectangle has empty interior. */
      bool empty_interior() const { return _position.empty_interior(); }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;

      friend bool interiors_intersect<> (const GridRectangle<R>&, const GridRectangle<R>&);
      friend bool interiors_intersect<> (const GridRectangle<R>&, const GridMaskSet<R>&);
      friend bool subset<> (const GridRectangle<R>&, const GridMaskSet<R>&);
     private:
      const Grid<R>& _grid;
      LatticeRectangle _position;
    };

    template<class Base, class Value>
    class GridSetIterator 
      : public boost::iterator_adaptor<GridSetIterator<Base,Value>,Base,Value,boost::use_default,Value>
    { 
     public:
      typedef typename Value::real_type real_type;
      GridSetIterator(const Grid<real_type>& g, Base i) 
        : GridSetIterator::iterator_adaptor_(i), _grid(g) { }
     private:
      friend class boost::iterator_core_access;
      Value dereference() const { return Value(_grid,*this->base_reference()); }
      const Grid<real_type>& _grid;
    };


    template<class R>
    class GridMaskSetIterator {
      typedef GridMaskSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridCell<R> value_type;
      typedef GridCell<R> reference;
      typedef const GridCell<R>* pointer;
      typedef int difference_type;
     public:
      GridMaskSetIterator(const Grid<R>& g, LatticeMaskSetIterator iter);
      GridCell<R> operator*() const { return this->dereference(); }
      Self& operator++() { this->increment(); return *this; }
      Self operator++(int) { Self tmp=*this; this->increment(); return tmp; }
      bool operator==(const Self& other) const { return this->equal(other); }
      bool operator!=(const Self& other) const { return !this->equal(other); }
     private:
      bool equal(const GridMaskSetIterator& other) const { return this->_iter==other._iter; }
      void increment() { ++_iter; }
      GridCell<R> dereference() const { return GridCell<R>(*_grid_ptr,*_iter); }
     private:
      const Grid<R>* _grid_ptr;
      LatticeMaskSet::const_iterator _iter;
    };

    
    
    /*! \brief A denotable set on a finite grid, defined using a mask. */
    template<typename R>
    class GridMaskSet {
      friend class GridCellListSet<R>;
      friend class PartitionTreeSet<R>;
     public:
      typedef R real_type;
      typedef Point<R> state_type;
      typedef GridCell<R> value_type;
      typedef GridSetIterator< LatticeMaskSet::const_iterator, GridCell<R> > iterator;
      typedef GridSetIterator< LatticeMaskSet::const_iterator, GridCell<R> > const_iterator;

      /*!\brief Construct an empty set from a finite grid. */
      GridMaskSet(const FiniteGrid<R>& g);
     
      /*!\brief Construct an empty set from a grid, and a bounding box. */
      GridMaskSet(const Grid<R>& g, const LatticeRectangle& b);
     
      /*!\brief Construct a set from a grid, and a lattice mask set. */
      GridMaskSet(const Grid<R>& g, const LatticeMaskSet& ms);
     
      /*!\brief Construct a set from a grid, a bounding box, and a mask */
      GridMaskSet(const Grid<R>& g, const LatticeRectangle& b, const BooleanArray& m);

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
        if(this->grid()==gms.grid() && this->bounds()==gms.bounds()) {
          return this->mask()==gms.mask();
        }
        throw(std::domain_error("Can only compare GridMaskSets on the same grid"));
      }

      /*! \brief Inequality operator. */
      bool operator!=(const GridMaskSet<R>& gms) { return !(*this==gms); }

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return *_grid_ptr; }

      /*! \brief The underlying grid. */
      const LatticeMaskSet& lattice_set() const { return _lattice_set; }

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const { return _lattice_set.dimension(); }

       /*! \brief The number of elements in the mask. */
      size_type capacity() const { return _lattice_set.capacity(); }
      
      /*! \brief The lowest position in the grid. */
      const IndexArray& lower() const { return _lattice_set.lower(); }

      /*! \brief The highest position in the grid. */
      const IndexArray& upper() const { return _lattice_set.upper(); }

      /*! \brief The bounding rectangle in the grid. */
      const LatticeRectangle& bounds() const { return _lattice_set.bounds(); }

      /*! \brief The bounding rectangle in the grid. */
      const SizeArray& sizes() const { return _lattice_set.sizes(); }

      /*! \brief The bounding rectangle in the grid. */
      const SizeArray& strides() const { return _lattice_set.strides(); }

      /*! \brief The number of cells in the grid. */
      const BooleanArray& mask() const { return _lattice_set.mask(); }

      /*! \brief Returns true if the set is empty. */
      bool empty() const { return _lattice_set.empty(); }
      
      /*! \brief The number of cells in the grid. */
      size_type size() const { return _lattice_set.size(); }
      
      /*! \brief The ith nonempty cell in the grid. */
      GridCell<R> operator[](size_type i) const { return GridCell<R>(*_grid_ptr,_lattice_set[i]); }

      /*! \brief A constant iterator to the beginning of the set. */
      const_iterator begin() const { return const_iterator(*this->_grid_ptr,this->_lattice_set.begin()); }
      /*! \brief A constant iterator to the end of the set. */
      const_iterator end() const { return const_iterator(*this->_grid_ptr,this->_lattice_set.end()); }

      /*! \brief The rectangle bounding the region of the mask. */
      GridRectangle<R> bounding_box() const { return GridRectangle<R>(grid(),bounds()); }
      
      /*! \brief Removes all cells. */
      void clear() {
        _lattice_set.clear();
      }

      /*! \brief Adjoins a cell to the set. */
      void adjoin(const GridCell<R>& c) {
        assert(c.grid()==this->grid());
        this->_lattice_set.adjoin(c._position);
      }

      /*! \brief Adjoins a rectangle to the set. */
      void adjoin(const GridRectangle<R>& r) {
        assert(r.grid()==this->grid());
        this->_lattice_set.adjoin(r._position);
      }

      /*! \brief Adjoins a GridMaskSet to the set. */
      void adjoin(const GridMaskSet<R>& gms) {
        //std::cerr << "GridMaskSet<R>::adjoin(const GridMaskSet<R>&)" << std::endl;
        assert(gms.grid()==this->grid());
        this->_lattice_set.adjoin(gms._lattice_set);
      }

      /*! \brief Adjoins a GridCellListSet to the set. */
      void adjoin(const GridCellListSet<R>& cls) {
        assert(cls.grid()==this->grid());
        this->_lattice_set.adjoin(cls._lattice_set);
      }

      /*! \brief Adjoins a GridRectangleListSet to the set. */
      void adjoin(const GridRectangleListSet<R>& rls) {
        assert(rls.grid()==this->grid());
        this->_lattice_set.adjoin(rls._lattice_set);
      }

      GridMaskSet neighbourhood() const {
        return GridMaskSet(this->grid(),this->_lattice_set.neighbourhood());
      }

      GridMaskSet adjoining() const {
        return GridMaskSet(this->grid(),this->_lattice_set.adjoining());
      }

      friend bool interiors_intersect<> (const GridRectangle<R>&, const GridMaskSet<R>&);
      friend bool interiors_intersect<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend bool subset<> (const GridRectangle<R>&, const GridMaskSet<R>&);
      friend bool subset<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> join<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> regular_intersection<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> difference<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      LatticeMaskSet _lattice_set;
    };



    template<class R>
    class GridCellListSetIterator {
      typedef GridCellListSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridCell<R> value_type;
      typedef GridCell<R> reference;
      typedef const GridCell<R>* pointer;
      typedef int difference_type;
     public:
      GridCellListSetIterator(const Grid<R>& g, LatticeCellListSetIterator iter);
      GridCell<R> operator*() const { return this->dereference(); }
      Self& operator++() { this->increment(); return *this; }
      Self operator++(int) { Self tmp=*this; this->increment(); return tmp; }
      bool operator==(const Self& other) const { return this->equal(other); }
      bool operator!=(const Self& other) const { return !this->equal(other); }
     private:
      bool equal(const GridCellListSetIterator& other) const { return this->_iter==other._iter; }
      void increment() { ++_iter; }
      GridCell<R> dereference() const { return GridCell<R>(*_grid_ptr,*_iter); }
     private:
      const Grid<R>* _grid_ptr;
      LatticeCellListSet::const_iterator _iter;
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
      typedef GridCell<R> value_type;
      typedef GridSetIterator< LatticeCellListSet::const_iterator, GridCell<R> > iterator;
      typedef GridSetIterator< LatticeCellListSet::const_iterator, GridCell<R> > const_iterator;

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
      dimension_type dimension() const { return _lattice_set.dimension(); }

      /*! \brief True if the set is empty. */
      bool empty() const { return _lattice_set.empty(); }

      /*! \brief The numeber of cells in the list. */
      size_type size() const { return _lattice_set.size(); }

      /*! \brief The lattice set. */
      const LatticeCellListSet& lattice_set() const { return _lattice_set; }

      /*!\brief The @a i th cell in the list. */
      GridCell<R> operator[] (const size_type i) const { return GridCell<R>(grid(),_lattice_set[i]); }

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const { return const_iterator(*_grid_ptr,_lattice_set.begin()); }

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const { return const_iterator(*_grid_ptr,_lattice_set.end()); }

      /*!\brief Append a GridCell to the list. */
      void adjoin(const GridCell<R>& c) { _lattice_set.adjoin(c.position()); }

      friend std::ostream& operator<< <> (std::ostream&, const GridCellListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      LatticeCellListSet _lattice_set;
    };





    template<class R>
    class GridRectangleListSetIterator {
      typedef GridRectangleListSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridRectangle<R> value_type;
      typedef GridRectangle<R> reference;
      typedef const GridRectangle<R>* pointer;
      typedef int difference_type;
     public:
      GridRectangleListSetIterator(const Grid<R>& g, LatticeRectangleListSetIterator iter);
      GridRectangle<R> operator*() const { return this->dereference(); }
      Self& operator++() { this->increment(); return *this; }
      Self operator++(int) { Self tmp=*this; this->increment(); return tmp; }
      bool operator==(const Self& other) const { return this->equal(other); }
      bool operator!=(const Self& other) const { return !this->equal(other); }
     private:
      bool equal(const GridRectangleListSetIterator& other) const { return this->_iter==other._iter; }
      void increment() { ++_iter; }
      GridRectangle<R> dereference() const { return GridRectangle<R>(*_grid_ptr,*_iter); }
     private:
      const Grid<R>* _grid_ptr;
      LatticeRectangleListSet::const_iterator _iter;
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
      typedef GridRectangle<R> value_type;
      typedef GridSetIterator< LatticeRectangleListSet::const_iterator, GridRectangle<R> > iterator;
      typedef GridSetIterator< LatticeRectangleListSet::const_iterator, GridRectangle<R> > const_iterator;

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
      dimension_type dimension() const { return _lattice_set.dimension(); }

      /*! \brief True if the set is empty. */
      bool empty() const { return _lattice_set.empty(); }

      /*! \brief The number of rectangles in the list. */
      size_type size() const { return _lattice_set.size(); }

      /*! \brief The rectangles on the underlying lattice. */
      const LatticeRectangleListSet& lattice_set() { return _lattice_set; }
      
      /*!\brief Return the @a i th rectangle in the list. */
      GridRectangle<R> operator[] (const size_t i) const {
        return GridRectangle<R>(grid(),_lattice_set[i]);
      }

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const { return const_iterator(*_grid_ptr,_lattice_set.begin()); }

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const { return const_iterator(*_grid_ptr,_lattice_set.end()); }

      /*!\brief Append a GridRectangle to the list. */
      void adjoin(const GridRectangle<R>& r) {
        _lattice_set.adjoin(r._position);
      }

      friend std::ostream& operator<< <> (std::ostream&, const GridRectangleListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      LatticeRectangleListSet _lattice_set;
    };


  }
}

#endif /* _GRID_SET_H */
