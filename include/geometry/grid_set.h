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

#include <boost/iterator/iterator_adaptor.hpp>

#include "../declarations.h"

#include "../base/array.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/grid.h"

/*TODO: Unify bounds in FiniteGrid, and make GridMaskSet use only finite grid bounds*/

namespace Ariadne {
  namespace Geometry {

    template<typename R> bool interiors_intersect(const GridBlock<R>&, const GridBlock<R>&);
    template<typename R> bool interiors_intersect(const GridBlock<R>&, const GridMaskSet<R>&);
    template<typename R> bool interiors_intersect(const GridMaskSet<R>&, const GridBlock<R>&);
    template<typename R> bool interiors_intersect(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool interiors_intersect(const Rectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool interiors_intersect(const GridMaskSet<R>& gm, const Rectangle<R>& r);

    template<typename R> bool disjoint(const GridBlock<R>&, const GridBlock<R>&);
    template<typename R> bool disjoint(const GridBlock<R>&, const GridMaskSet<R>&);
    template<typename R> bool disjoint(const GridMaskSet<R>&, const GridBlock<R>&);
    template<typename R> bool disjoint(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool disjoint(const Rectangle<R>&, const GridMaskSet<R>&);
    template<typename R> bool disjoint(const GridMaskSet<R>& gm, const Rectangle<R>& r);

    template<typename R> bool subset(const GridBlock<R>&, const GridBlock<R>&);
    template<typename R> bool subset(const GridCell<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridBlock<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridCellListSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridBlockListSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool subset(const Rectangle<R>&, const GridBlock<R>&);
    template<typename R> bool subset(const Rectangle<R>&, const GridMaskSet<R>&);

    template<typename R> bool inner_subset(const GridBlock<R>&, const GridBlock<R>&);
    template<typename R> bool inner_subset(const GridBlock<R>&, const GridMaskSet<R>&);
    template<typename R> bool inner_subset(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> bool inner_subset(const Rectangle<R>&, const GridBlock<R>&);
    template<typename R> bool inner_subset(const Rectangle<R>&, const GridMaskSet<R>&);

    template<typename R> GridMaskSet<R> regular_intersection(const GridBlock<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> regular_intersection(const GridMaskSet<R>&, const GridBlock<R>&);
    template<typename R> GridCellListSet<R> regular_intersection(const GridCellListSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridCellListSet<R> regular_intersection(const GridMaskSet<R>&, const GridCellListSet<R>&);
    template<typename R> GridCellListSet<R> difference(const GridCellListSet<R>&, const GridMaskSet<R>&);

    template<typename R> GridMaskSet<R> regular_intersection(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> join(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<typename R> GridMaskSet<R> difference(const GridMaskSet<R>&, const GridMaskSet<R>&);

    template<typename R>
    GridBlock<R>
    outer_approximation(const Rectangle<R>& p, const Grid<R>& g);

    template<typename R>
    GridBlock<R>
    inner_approximation(const Rectangle<R>& p, const Grid<R>& g);


    template<typename R>
    GridBlock<R>
    over_approximation(const Rectangle<R>& p, const Grid<R>& g);

    template<typename R>
    GridCellListSet<R>
    over_approximation(const Zonotope<R>& p, const Grid<R>& g);

    template<typename R>
    GridCellListSet<R>
    over_approximation(const Polytope<R>& p, const Grid<R>& g);
    
    
    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g); 

    template<typename R>
    GridMaskSet<R>
    over_approximation(const GridMaskSet<R>& gm, const FiniteGrid<R>& g);
    
    template<typename R>
    GridMaskSet<R>
    over_approximation(const PartitionTreeSet<R>& pts, const FiniteGrid<R>& g);
    
    
    
    template<typename R>
    GridBlock<R>
    under_approximation(const Rectangle<R>& p, const Grid<R>& g);

    template<typename R>
    GridCellListSet<R>
    under_approximation(const Zonotope<R>& p, const Grid<R>& g);
   
    template<typename R>
    GridCellListSet<R>
    under_approximation(const Polytope<R>& p, const Grid<R>& g);
       
    
    template<typename R>
    GridMaskSet<R>
    under_approximation(const ListSet<R,Rectangle>& ls, const FiniteGrid<R>& g); 

    template<typename R>
    GridMaskSet<R>
    under_approximation(const GridMaskSet<R>& gm, const FiniteGrid<R>& g);
    
    template<typename R>
    GridMaskSet<R>
    under_approximation(const PartitionTreeSet<R>& gm, const FiniteGrid<R>& g);
    
    

    template<typename R> std::ostream& operator<<(std::ostream&, const GridCell<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridBlock<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridBlockListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridCellListSet<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const GridMaskSet<R>&);
    
    /*! \brief A cell in a grid.
     *  \ingroup BasicSet
     *  \ingroup Grid
     */
    template<typename R>
    class GridCell 
      : public RectangleExpression< GridCell<R> >
    {
      friend class GridBlock<R>;
      friend class GridMaskSet<R>;
      friend class GridCellListSet<R>;
     public:
      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      
      /*!\brief Construct from a grid and an unit grid cell. */
      GridCell(const Grid<R>& g, const Combinatoric::LatticeCell& pos);

      /*!\brief Construct from a grid and an unit grid cell. */
      GridCell(const Grid<R>& g, const IndexArray& pos);

      /*!\brief The grid containing the cell. */
      const Grid<R>& grid() const { return _grid; }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return _lattice_set.dimension(); }

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;

      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;

      /*!\brief The position of the cell in the grid. */
      const Combinatoric::LatticeCell& lattice_set() const { return _lattice_set; }

      /*!\brief A rectangle containing the grid rectangle. */
      Rectangle<R> bounding_box() const { return *this; }
     private:
      const Grid<R>& _grid;
      Combinatoric::LatticeCell _lattice_set;
    };


    /*! \brief A rectangle in a grid.
     *  \ingroup BasicSet
     *  \ingroup Grid
     */
    template<typename R>
    class GridBlock 
      : public RectangleExpression< GridBlock<R> >
    {
      friend class GridCell<R>;
      friend class GridMaskSet<R>;
      friend class GridBlockListSet<R>;
     public:
      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;

      /*!\brief Construct an empty rectangle on a grid. */
      GridBlock(const Grid<R>& g);
      /*!\brief Construct from a grid and a bounding block. */
      GridBlock(const Grid<R>& g, const Combinatoric::LatticeBlock& b);
      /*!\brief Construct from a grid and two integer arrays giving the corners. */
      GridBlock(const Grid<R>& g, const IndexArray& l, const IndexArray& u);
      /*!\brief Construct from a grid and an ordinary rectangle. */
      GridBlock(const Grid<R>& g, const Rectangle<R>& r);
      /*!\brief Construct from a GridCell. */
      GridBlock(const GridCell<R>& c);

      /*!\brief The grid containing the rectangle. */
      const Grid<R>& grid() const { return _grid; }

      /*!\brief The dimension of the rectangle. */
      dimension_type dimension() const { return _lattice_set.dimension(); }

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;
      
      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;
      
      /*!\brief The position of the rectangle in the grid. */
      const Combinatoric::LatticeBlock& lattice_set() const { return _lattice_set; }

      /*!\brief Tests if the rectangle is empty. */
      bool empty() const { return _lattice_set.empty(); }
      /*!\brief Tests if the rectangle has empty interior. */
      bool empty_interior() const { return _lattice_set.empty_interior(); }

      /*!\brief A rectangle containing the grid rectangle. */
      Rectangle<R> bounding_box() const { return *this; }

      friend bool interiors_intersect<> (const GridBlock<R>&, const GridBlock<R>&);
      friend bool interiors_intersect<> (const GridBlock<R>&, const GridMaskSet<R>&);
      friend bool interiors_intersect<> (const GridMaskSet<R>&, const GridBlock<R>&);
      friend bool subset<> (const GridBlock<R>&, const GridMaskSet<R>&);
     private:
      const Grid<R>& _grid;
      Combinatoric::LatticeBlock _lattice_set;
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
    class GridCellListSetIterator {
      typedef GridCellListSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridCell<R> value_type;
      typedef GridCell<R> reference;
      typedef const GridCell<R>* pointer;
      typedef int difference_type;
     public:
      GridCellListSetIterator(const Grid<R>& g, Combinatoric::LatticeCellListSetIterator iter);
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
      Combinatoric::LatticeCellListSet::const_iterator _iter;
    };


    /*! \brief A denotable set on a grid, defined using a list of cells.
     *  \ingroup DenotableSet
     *  \ingroup Grid
     */
    template<typename R>
    class GridCellListSet {
      friend class GridBlockListSet<R>;
      friend class GridMaskSet<R>;
     public:
      typedef GridSetIterator< Combinatoric::LatticeCellListSet::const_iterator, GridCell<R> > iterator;
      typedef GridSetIterator< Combinatoric::LatticeCellListSet::const_iterator, GridCell<R> > const_iterator;

      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      /*! \brief The type of basic set contained by the denotable set. */
      typedef GridCell<R> value_type;
      
      /*!\brief Construct an empty set based on a Grid. */
      GridCellListSet(const Grid<R>& g);

      /*!\brief Construct from a grid and a lattice cell list set. */
      GridCellListSet(const Grid<R>& g, const Combinatoric::LatticeCellListSet& lcls);

      /*!\brief Construct from a GridMaskSet. */
      GridCellListSet(const GridMaskSet<R>& gms);

      /*!\brief Copy constructor. */
      GridCellListSet(const GridCellListSet<R>& gcls);

      /*!\brief Construct from a GridBlockListSet. */
      GridCellListSet(const GridBlockListSet<R>& grls);

      /*!\brief Construct from a GridBlockListSet. */
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
      const Combinatoric::LatticeCellListSet& lattice_set() const { return _lattice_set; }

      /*!\brief The @a i th cell in the list. */
      GridCell<R> operator[] (const size_type i) const { return GridCell<R>(grid(),_lattice_set[i]); }

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const { return const_iterator(*_grid_ptr,_lattice_set.begin()); }

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const { return const_iterator(*_grid_ptr,_lattice_set.end()); }

      /*!\brief Append a GridCell to the list. */
      void adjoin(const GridCell<R>& c) { _lattice_set.adjoin(c.lattice_set()); }
      /*!\brief Append a GridBlock to the list. */
      void adjoin(const GridBlock<R>& bl) { _lattice_set.adjoin(bl.lattice_set()); }
      /*!\brief Append a GridCellListSet to the list. */
      void adjoin(const GridCellListSet<R>& cls) { _lattice_set.adjoin(cls.lattice_set()); }

      /*! \brief Empties the set. */
      void clear();

      friend std::ostream& operator<< <> (std::ostream&, const GridCellListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      Combinatoric::LatticeCellListSet _lattice_set;
    };





    template<class R>
    class GridBlockListSetIterator {
      typedef GridBlockListSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridBlock<R> value_type;
      typedef GridBlock<R> reference;
      typedef const GridBlock<R>* pointer;
      typedef int difference_type;
     public:
      GridBlockListSetIterator(const Grid<R>& g, Combinatoric::LatticeBlockListSetIterator iter);
      GridBlock<R> operator*() const { return this->dereference(); }
      Self& operator++() { this->increment(); return *this; }
      Self operator++(int) { Self tmp=*this; this->increment(); return tmp; }
      bool operator==(const Self& other) const { return this->equal(other); }
      bool operator!=(const Self& other) const { return !this->equal(other); }
     private:
      bool equal(const GridBlockListSetIterator& other) const { return this->_iter==other._iter; }
      void increment() { ++_iter; }
      GridBlock<R> dereference() const { return GridBlock<R>(*_grid_ptr,*_iter); }
     private:
      const Grid<R>* _grid_ptr;
      Combinatoric::LatticeBlockListSet::const_iterator _iter;
    };


    /*! \brief A denotable set on a grid, defined using a list of rectangles.
     *  \ingroup DenotableSet
     *  \ingroup Grid
     */
    template<typename R>
    class GridBlockListSet {
      friend class GridMaskSet<R>;
      friend class GridCellListSet<R>;
     public:
      typedef GridSetIterator< Combinatoric::LatticeBlockListSet::const_iterator, GridBlock<R> > iterator;
      typedef GridSetIterator< Combinatoric::LatticeBlockListSet::const_iterator, GridBlock<R> > const_iterator;

      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      /*! \brief The type of basic set contained by the denotable set. */
      typedef GridBlock<R> value_type;
      /*!\brief Destructor. */
      ~GridBlockListSet() { 
        //delete _grid_ptr; }
      }
      
      /*!\brief Construct from a Grid. */
      GridBlockListSet(const Grid<R>& g);

      /*!\brief Construct from a Grid and a lattice block list set. */
      GridBlockListSet(const Grid<R>& g, const Combinatoric::LatticeBlockListSet& lbls);

      /*!\brief Construct from a ListSet of Rectangles. */
      explicit GridBlockListSet(const ListSet<R,Rectangle>& s);

      /*!\brief Construct from a PartitionTreeSet. */
      /* FIXME: This constructor is only included since Boost Python doesn't find the conversion operator */
      explicit GridBlockListSet(const PartitionTreeSet<R>& s);

      /*!\brief Copy constructor. */
      GridBlockListSet(const GridBlockListSet<R>& s);

      /*!\brief Construct a set on a finer grid. */
      GridBlockListSet(const Grid<R>& g, const ListSet<R,Rectangle>& s);

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
      const Combinatoric::LatticeBlockListSet& lattice_set() const { return _lattice_set; }
      
      /*!\brief Return the @a i th rectangle in the list. */
      GridBlock<R> operator[] (const size_t i) const {
        return GridBlock<R>(grid(),_lattice_set[i]);
      }

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const { return const_iterator(*_grid_ptr,_lattice_set.begin()); }

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const { return const_iterator(*_grid_ptr,_lattice_set.end()); }

      /*!\brief Append a GridBlock to the list. */
      void adjoin(const GridBlock<R>& r) {
        _lattice_set.adjoin(r._lattice_set);
      }
      
      /*! \brief Empties the set. */
      void clear();

      friend std::ostream& operator<< <> (std::ostream&, const GridBlockListSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      Combinatoric::LatticeBlockListSet _lattice_set;
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
      GridMaskSetIterator(const Grid<R>& g, Combinatoric::LatticeMaskSetIterator iter);
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
      Combinatoric::LatticeMaskSet::const_iterator _iter;
    };

    
    
    /*! \brief A denotable set on a finite grid, defined using a mask. 
     *  \ingroup DenotableSet
     *  \ingroup Grid
     */
    template<typename R>
    class GridMaskSet {
      friend class GridCellListSet<R>;
      friend class PartitionTreeSet<R>;
     public:
      typedef GridSetIterator< Combinatoric::LatticeMaskSet::const_iterator, GridCell<R> > iterator;
      typedef GridSetIterator< Combinatoric::LatticeMaskSet::const_iterator, GridCell<R> > const_iterator;

      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      /*! \brief The type of basic set contained by the denotable set. */
      typedef GridCell<R> value_type;
      /*! \brief The type of object returned by indexing. */

      /*!\brief Construct an empty set from a finite grid. (Deprecated)
       *
       * \deprecated Use GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b) instead.
       */
      GridMaskSet(const FiniteGrid<R>& g);
     
      /*!\brief Construct a set from a finite grid and a mask. (Deprecated)
       *
       * \deprecated Use GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b, const BooleanArray& m) instead.
       */
      GridMaskSet(const FiniteGrid<R>& g, const BooleanArray& m);

      /*!\brief Construct an empty set based on grid \a g and with cells in the block \a b. */
      GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b);
     
      /*!\brief Construct a set from a grid, a lattice rectangle and a mask. */
      GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b, const BooleanArray& m);
     
      /*!\brief Construct a set from a grid, and a lattice mask set. */
      GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeMaskSet& ms);
     
      /*!\brief Copy constructor. */
      GridMaskSet(const GridMaskSet<R>& gms);

      /*!\brief Construct from a %ListSet of %Rectangle. */
      GridMaskSet(const ListSet<R,Rectangle>& ls);

      /*!\brief Construct from a %GridCellListSet. */
      GridMaskSet(const GridCellListSet<R>& gcls);

      /*!\brief Construct from a %GridBlockListSet. */
      GridMaskSet(const GridBlockListSet<R>& grls);


      /*!\brief Convert to a %ListSet of BS. */
      template <template<typename> class BS> operator ListSet<R,BS> () const;
      operator ListSet<R,Rectangle> () const;
      
      /*!\brief Equality operator. (Deprecated)
       *
       * \deprecated Equality operator for sets is deprecated.
       */
      bool operator==(const GridMaskSet<R>& gms) {
        throw(std::domain_error("Equality operator for sets is deprecated"));
      }

      /*! \brief Inequality operator. (Deprecated)
       *
       * \deprecated Equality operator for sets is deprecated.
       */
      bool operator!=(const GridMaskSet<R>& gms) { return !(*this==gms); }

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const { return *_grid_ptr; }

      /*! \brief The bounds on elements in the set. */
      GridBlock<R> bounds() const { return GridBlock<R>(*_grid_ptr,_lattice_set.block()); }

      /*! \brief The underlying lattice set. */
      const Combinatoric::LatticeMaskSet& lattice_set() const { return _lattice_set; }

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const { return _lattice_set.dimension(); }

       /*! \brief The number of elements in the mask. */
      size_type capacity() const { return _lattice_set.capacity(); }
      
      /*! \brief The block of cells available in the lattice. */
      const Combinatoric::LatticeBlock& block() const { return _lattice_set.block(); }

      /*! \brief The lowest position in the grid. (Deprecated)
       *
       * \deprecated Use block().lower() instead.
       */
      //const IndexArray& lower() const { return _lattice_set.block().lower(); }

      /*! \brief The highest position in the grid. (Deprecated)
       *
       * \deprecated Use block().lower() instead.
       */
      //const IndexArray& upper() const { return _lattice_set.block().upper(); }

      /*! \brief The number of cells in each dimension. */
      const SizeArray& sizes() const { return _lattice_set.sizes(); }

      /*! \brief The number of c. 
       *
       * \deprecated Use block().sizes() instead.
       */
      //const SizeArray& strides() const { return _lattice_set.block().strides(); }

      /*! \brief The mask giving the cells in the set. */
      const BooleanArray& mask() const { return _lattice_set.mask(); }

      /*! \brief Returns true if the set is empty. */
      bool empty() const { return _lattice_set.empty(); }
      
      /*! \brief Returns true if the set is empty. */
      bool bounded() const { return _lattice_set.bounded(); }
      
      /*! \brief The number of cells in the grid. */
      size_type size() const { return _lattice_set.size(); }
      
      /*! \brief The ith nonempty cell in the grid. */
      GridCell<R> operator[](size_type i) const { return GridCell<R>(*_grid_ptr,_lattice_set[i]); }

      /*! \brief A constant iterator to the beginning of the set. */
      const_iterator begin() const { return const_iterator(*this->_grid_ptr,this->_lattice_set.begin()); }
      /*! \brief A constant iterator to the end of the set. */
      const_iterator end() const { return const_iterator(*this->_grid_ptr,this->_lattice_set.end()); }

      /*! \brief The rectangle bounding the region of the mask. */
      Rectangle<R> bounding_box() const { return GridBlock<R>(grid(),bounds()); }
     
      /*! \brief Removes all cells. */
      void clear();

      /*! \brief Adjoins a cell to the set. */
      void adjoin(const GridCell<R>& c) {
        assert(c.grid()==this->grid());
        this->_lattice_set.adjoin(c._lattice_set);
      }

      /*! \brief Adjoins a rectangle to the set. */
      void adjoin(const GridBlock<R>& r) {
        assert(r.grid()==this->grid());
        
        if (r.empty()) return;
        this->_lattice_set.adjoin(r._lattice_set);
      }

      /*! \brief Adjoins a GridMaskSet to the set. */
      void adjoin(const GridMaskSet<R>& gms) {
        assert(gms.grid()==this->grid());

        if (gms.empty()) 
          return;
        
        this->_lattice_set.adjoin(gms._lattice_set);
      }

      /*! \brief Adjoins a GridCellListSet to the set. */
      void adjoin(const GridCellListSet<R>& cls) {
        assert(cls.grid()==this->grid());
        
        if (cls.empty()) 
          return;
        
        this->_lattice_set.adjoin(cls._lattice_set);
      }

      /*! \brief Adjoins a GridBlockListSet to the set. */
      void adjoin(const GridBlockListSet<R>& rls) {
        assert(rls.grid()==this->grid());
        
        if (rls.empty()) return;
        this->_lattice_set.adjoin(rls._lattice_set);
      }

      /*!\brief The one-box neighbourhood, consisting of all cells whose closure intersects the set. */
      GridMaskSet neighbourhood() const {
        return GridMaskSet(this->grid(),this->_lattice_set.neighbourhood());
      }

      /*!\brief The set of all cells which share a facet a cell in the set. */
      GridMaskSet adjoining() const {
        return GridMaskSet(this->grid(),this->_lattice_set.adjoining());
      }

      friend bool interiors_intersect<> (const GridBlock<R>&, const GridMaskSet<R>&);
      friend bool interiors_intersect<> (const GridMaskSet<R>&, const GridBlock<R>&);
      friend bool interiors_intersect<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend bool subset<> (const GridBlock<R>&, const GridMaskSet<R>&);
      friend bool subset<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> join<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> regular_intersection<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> difference<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
     private:
      const Grid<R>* _grid_ptr;
      Combinatoric::LatticeMaskSet _lattice_set;
    };

    template<typename R>
    template <template<typename> class BS>
    inline
    GridMaskSet<R>::operator ListSet<R,BS> () const 
    {
      ListSet<R,BS> result(this->dimension());
      Rectangle<R> r(this->dimension());
      BS<R> bs(this->dimension());
      for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        r=iter;
        bs=BS<R>(r);
        result.push_back(bs);
      }
      return result;
    }

  }
}

#endif /* _GRID_SET_H */
