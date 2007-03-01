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

#ifndef ARIADNE_GRID_SET_H
#define ARIADNE_GRID_SET_H

#include <iosfwd>

#include <boost/iterator/iterator_adaptor.hpp>

#include "../declarations.h"

#include "../base/array.h"
#include "../base/tribool.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/grid.h"

/*TODO: Unify bounds in FiniteGrid, and make GridMaskSet use only finite grid bounds*/

namespace Ariadne {
  namespace Geometry {

    template<class Base, class Value> class GridSetIterator;
    template<class R> class GridBlockIterator;
    template<class R> class GridCellListSetIterator;
    template<class R> class GridMaskSetIterator;

    /*! \brief A unit cell in a grid.
     *
     *  A %GridCell is defined by mapping a LatticeCell \a lc into \f$\mathbb{R}^d\f$
     *  via a grid \a g. The bounds of the ith coordinate of the cell are given by
     *  g.subdivision_coordinate(i,lc.lower_bound(i)) and
     *  g.subdivision_coordinate(i,lc.upper_bound(i)).
     *
     *  A %GridCell satisfies the requirements of a RectangleExpression.
     *  \ingroup BasicSet
     *  \ingroup Grid
     */
    template<class R>
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
      const Grid<R>& grid() const;

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const;

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;

      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;

      /*!\brief The position of the cell in the grid. */
      const Combinatoric::LatticeCell& lattice_set() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*!\brief A rectangle containing the grid cell. */
      Rectangle<R> bounding_box() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private: 
      static void _instantiate_geometry_operators();
     private:
      const Grid<R>& _grid_ref;
      Combinatoric::LatticeCell _lattice_set;
    };


    /*! \brief A rectangular block in a grid.
     *
     *  A %GridBlock is defined by mapping a LatticeBlock \a lb into \f$\mathbb{R}^d\f$
     *  via a Grid \a g. The bounds of the ith coordinate of the cell are given by
     *  g.subdivision_coordinate(i,lb.lower_bound(i)) and
     *  g.subdivision_coordinate(i,lb.upper_bound(i)).
     *
     *  A %GridBlock satisfies all the requirements of a RectangleExpression.
     *  \ingroup BasicSet
     *  \ingroup Grid
     */
    template<class R>
    class GridBlock 
      : public RectangleExpression< GridBlock<R> >
    {
      friend class GridCell<R>;
      friend class GridMaskSet<R>;
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
      GridBlock(const GridCell<R>& gc);
      
      /*!\brief Copy constructor. */
      GridBlock(const GridBlock<R>& gb);
      /*!\brief Copy assignment. */
      GridBlock<R>& operator=(const GridBlock<R>& gb);

      /*!\brief The grid containing the rectangle. */
      const Grid<R>& grid() const;

      /*!\brief The dimension of the rectangle. */
      dimension_type dimension() const;

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;
      
      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;
      
      /*!\brief The position of the rectangle in the grid. */
      const Combinatoric::LatticeBlock& lattice_set() const;

      /*!\brief Tests if the rectangle is empty. */
      tribool empty() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*!\brief A rectangle containing the grid rectangle. */
      Rectangle<R> bounding_box() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private: 
      static void _instantiate_geometry_operators();
     private:
      const Grid<R>& _grid_ref;
      Combinatoric::LatticeBlock _lattice_set;
    };



    /*! \brief A denotable set on a grid, defined as a list of cells.
     *
     *  The data for a %GridCellListSet provided by a Grid and a LatticeCellList.
     *
     *  A cell list set is useful for storing small sparse sets, such as an 
     *  over-approximation to a basic set.
     *
     *  \ingroup DenotableSet
     *  \ingroup Grid
     */
    template<class R>
    class GridCellListSet 
    {
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

      /*!\brief Copy assignment. */
      GridCellListSet<R>& operator=(const GridCellListSet<R>& gcls);

      /*!\brief Convert to a ListSet of Rectangles. */
      operator ListSet< Rectangle<R> >() const;

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const;

      /*! \brief The space dimension of the set. */
      dimension_type dimension() const;

      /*! \brief True if the set is empty. */
      tribool empty() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*! \brief The numeber of cells in the list. */
      size_type size() const;

      /*! \brief Empties the set. */
      void clear();

      /*! \brief The lattice set. */
      const Combinatoric::LatticeCellListSet& lattice_set() const;

      /*!\brief The @a i th cell in the list. */
      GridCell<R> operator[] (const size_type i) const;

      /*! \brief A constant iterator to the beginning of the list. */
      const_iterator begin() const;

      /*! \brief A constant iterator to the end of the list. */
      const_iterator end() const;

      /*! \brief Sorts the cells lexicographically, removing duplicates. */
      void unique_sort(); 

      /*!\brief Append a GridCell to the list. */
      void adjoin(const GridCell<R>& c);
      /*!\brief Append a GridBlock to the list. */
      void adjoin(const GridBlock<R>& bl);
      /*!\brief Append a GridCellListSet to the list. */
      void adjoin(const GridCellListSet<R>& cls);

      /*! \brief Adjoins an over-approximation of the set \a s. */
      template<class S>
      void adjoin_over_approximation(const S& s);

      /*! \brief Adjoins an under-approximation of the set \a s. */
      template<class S>
      void adjoin_under_approximation(const S& s);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private: 
      static void _instantiate_geometry_operators();
     private:
      const Grid<R>& _grid_ref;
      Combinatoric::LatticeCellListSet _lattice_set;
    };



    
    
    /*! \brief A denotable set on a finite grid, defined using a mask over a block of cells.
     *
     *  A %GridMaskSet is defined using a Grid and a LatticeMaskSet.
     * 
     *  A %GridMaskSet is especially useful for storing large sets with a lot of
     *  fine structure. Elementary set operations such as intersection, union and
     *  set difference can be implemented highly efficiently for GridMaskSets 
     *  defined over the same block of cells.
     * 
     *  For structured sets, a PartitionTreeSet is usually more efficient.
     *
     *  \ingroup DenotableSet
     *  \ingroup Grid
     */
    template<class R>
    class GridMaskSet 
      : public SetInterface<R>
    {
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
     
      /*!\brief Construct from a %GridCellListSet. */
      GridMaskSet(const GridCellListSet<R>& gcls);

      /*!\brief Copy constructor. */
      GridMaskSet(const GridMaskSet<R>& gms);

      /*!\brief Copy assignment. */
      GridMaskSet<R>& operator=(const GridMaskSet<R>& gms);

      /*!\brief Construct from a list set of rectangles. */
      explicit GridMaskSet(const ListSet< Rectangle<R> >& rls);

      //@{
      //! \name SetInterface methods
      /*! \brief Tests for disjointness with a Rectangle. */
      virtual GridMaskSet<R>* clone() const;

      /*! \brief The space dimension of the set. */
      virtual dimension_type dimension() const;

      /*!\brief Checks if a denotable set includes a point. */
      virtual tribool contains(const Point<R>& p) const;

      /*! \brief Tests for disjointness with a Rectangle. */
      virtual tribool disjoint(const Rectangle<R>& r) const;

      /*! \brief Tests for superset of a Rectangle. */ 
      virtual tribool superset(const Rectangle<R>& r) const;

      /*! \brief Tests for subset of a Rectangle. */
      virtual tribool subset(const Rectangle<R>& r) const;

      /*! \brief The rectangle bounding the region of the mask. */
      virtual Rectangle<R> bounding_box() const;
      //@}

      /*!\brief Convert to a %ListSet of BS. */
      template<class BS> operator ListSet<BS> () const;
      operator ListSet< Rectangle<R> > () const;
      
      /*!\brief Equality operator. (Deprecated)
       *
       * \deprecated Equality operator for sets is deprecated.
       */
      bool operator==(const GridMaskSet<R>& gms);

      /*! \brief Inequality operator. (Deprecated)
       *
       * \deprecated Equality operator for sets is deprecated.
       */
      bool operator!=(const GridMaskSet<R>& gms);

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const;

      /*! \brief The bounds on elements in the set. */
      GridBlock<R> bounds() const;

      /*! \brief The underlying lattice set. */
      const Combinatoric::LatticeMaskSet& lattice_set() const;

       /*! \brief The number of elements in the mask. */
      size_type capacity() const;
      
      /*! \brief The block of cells available in the lattice. */
      const Combinatoric::LatticeBlock& block() const;

      /*! \brief The number of cells in each dimension. */
      //const SizeArray& sizes() const;

      /*! \brief The mask giving the cells in the set. */
      const BooleanArray& mask() const;

      /*! \brief Returns true if the set is empty. */
      tribool empty() const;
      
      /*! \brief Returns true if the set is empty. */
      tribool bounded() const;
      
      /*! \brief The number of cells in the grid. */
      size_type size() const;
      
      /*! \brief The ith nonempty cell in the grid. */
      GridCell<R> operator[](size_type i) const;

      /*! \brief A constant iterator to the beginning of the set. */
      const_iterator begin() const;
      /*! \brief A constant iterator to the end of the set. */
      const_iterator end() const;

      /*! \brief Removes all cells. */
      void clear();

      /*! \brief Removes a cell from the set. */
      void remove(const GridCell<R>& gc);

      /*! \brief Adjoins a cell to the set. */
      void adjoin(const GridCell<R>& gc);

      /*! \brief Adjoins a rectangle to the set. */
      void adjoin(const GridBlock<R>& gb);

      /*! \brief Adjoins a GridCellListSet to the set. */
      void adjoin(const GridCellListSet<R>& gcls);

      /*! \brief Adjoins a GridMaskSet to the set. */
      void adjoin(const GridMaskSet<R>& gms);

      /*! \brief Intersects the set with a GridCellListSet. */
      void restrict(const GridCellListSet<R>& gcls);

      /*! \brief Intersects the set with a GridMaskSet. */
      void restrict(const GridMaskSet<R>& gms);

      /*!\brief The one-box neighbourhood, consisting of all cells whose closure intersects the set. */
      GridMaskSet neighbourhood() const;

      /*!\brief The set of all cells which share a facet a cell in the set. */
      GridMaskSet adjoining() const;

      /*! \brief Adjoins an over-approximation of the set \a s. */
      template<class S>
      void adjoin_over_approximation(const S& s);

      /*! \brief Adjoins an under-approximation of the set \a s. */
      template<class S>
      void adjoin_under_approximation(const S& s);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;

#ifdef DOXYGEN
      friend tribool overlap<> (const GridBlock<R>&, const GridMaskSet<R>&);
      friend tribool overlap<> (const GridMaskSet<R>&, const GridBlock<R>&);
      friend tribool overlap<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend tribool subset<> (const GridCell<R>&, const GridMaskSet<R>&);
      friend tribool subset<> (const GridBlock<R>&, const GridMaskSet<R>&);
      friend tribool subset<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> join<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> regular_intersection<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
      friend GridMaskSet<R> difference<> (const GridMaskSet<R>&, const GridMaskSet<R>&);
#endif
     private: 
      static void _instantiate_geometry_operators();
     private:
      const Grid<R>& _grid_ref;
      Combinatoric::LatticeMaskSet _lattice_set;
    };
    
    
    template<class R> tribool disjoint(const Rectangle<R>&, const GridMaskSet<R>&);
    template<class R> tribool disjoint(const GridMaskSet<R>& gm, const Rectangle<R>& r);
    template<class R> tribool subset(const Rectangle<R>&, const GridBlock<R>&);
    template<class R> tribool subset(const Rectangle<R>&, const GridMaskSet<R>&);
    template<class R> tribool subset(const GridMaskSet<R>&, const Rectangle<R>&);
   
    template<class R> tribool overlap(const GridBlock<R>&, const GridBlock<R>&);
    template<class R> tribool overlap(const GridBlock<R>&, const GridMaskSet<R>&);
    template<class R> tribool overlap(const GridMaskSet<R>&, const GridBlock<R>&);
    template<class R> tribool overlap(const GridMaskSet<R>&, const GridMaskSet<R>&);
    
    template<class R> tribool subset(const GridCell<R>&, const GridBlock<R>&);
    template<class R> tribool subset(const GridCell<R>&, const GridMaskSet<R>&);
    template<class R> tribool subset(const GridBlock<R>&, const GridBlock<R>&);
    template<class R> tribool subset(const GridCellListSet<R>&, const GridBlock<R>&);
    template<class R> tribool subset(const GridMaskSet<R>&, const GridBlock<R>&);
    template<class R> tribool subset(const GridCell<R>&, const GridMaskSet<R>&);
    template<class R> tribool subset(const GridBlock<R>&, const GridMaskSet<R>&);
    template<class R> tribool subset(const GridCellListSet<R>&, const GridMaskSet<R>&);
    template<class R> tribool subset(const GridMaskSet<R>&, const GridMaskSet<R>&);
    
    
    template<class R> GridCellListSet<R> regular_intersection(const GridCellListSet<R>&, const GridMaskSet<R>&);
    template<class R> GridCellListSet<R> regular_intersection(const GridMaskSet<R>&, const GridCellListSet<R>&);
    template<class R> GridMaskSet<R> regular_intersection(const GridBlock<R>&, const GridMaskSet<R>&);
    template<class R> GridMaskSet<R> regular_intersection(const GridMaskSet<R>&, const GridBlock<R>&);
    template<class R> GridMaskSet<R> regular_intersection(const GridMaskSet<R>&, const GridMaskSet<R>&);
    
    template<class R> GridMaskSet<R> difference(const GridMaskSet<R>&, const GridMaskSet<R>&);
    template<class R> GridCellListSet<R> difference(const GridCellListSet<R>&, const GridMaskSet<R>&);

    template<class R> GridMaskSet<R> join(const GridMaskSet<R>&, const GridMaskSet<R>&);

    
    
    template<class R> GridBlock<R> outer_approximation(const Rectangle<R>& r, const Grid<R>& g);
    
    
    template<class R> GridBlock<R> inner_approximation(const Rectangle<R>& r, const Grid<R>& g);
    
    
    template<class R> GridBlock<R> over_approximation(const Point< Interval<R> >& r, const Grid<R>& g);
    template<class R> GridBlock<R> over_approximation(const Rectangle<R>& r, const Grid<R>& g);
    
    template<class R> GridCellListSet<R> over_approximation(const Zonotope<R>& z, const Grid<R>& g);
    template<class R> GridCellListSet<R> over_approximation(const Polytope<R>& p, const Grid<R>& g);
    template<class R> GridCellListSet<R> over_approximation(const Polyhedron<R>& p, const Grid<R>& g);
    
    template<class R> GridCellListSet<R> over_approximation(const Zonotope< Interval<R> >& z, const Grid<R>& g);
      
    template<class R, class BS> GridMaskSet<R> over_approximation(const ListSet<BS>& ls, const FiniteGrid<R>& fg); 
    template<class R> GridMaskSet<R> over_approximation(const GridMaskSet<R>& gms, const FiniteGrid<R>& fg);
    template<class R> GridMaskSet<R> over_approximation(const PartitionTreeSet<R>& pts, const FiniteGrid<R>& fg);
    template<class R> GridMaskSet<R> over_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg);
    
    template<class R> GridMaskSet<R> over_approximation(const SetInterface<R>& set, const Grid<R>& g);
    
    
    template<class R> GridBlock<R> under_approximation(const Rectangle<R>& p, const Grid<R>& g);
    
    template<class R> GridCellListSet<R> under_approximation(const Zonotope<R>& p, const Grid<R>& g);
    template<class R> GridCellListSet<R> under_approximation(const Polytope<R>& p, const Grid<R>& g);
    template<class R> GridCellListSet<R> under_approximation(const Polyhedron<R>& p, const Grid<R>& g);
    
    template<class R> GridMaskSet<R> under_approximation(const ListSet< Rectangle<R> >& ls, const FiniteGrid<R>& fg); 
    template<class R> GridMaskSet<R> under_approximation(const GridMaskSet<R>& gms, const FiniteGrid<R>& fg);
    template<class R> GridMaskSet<R> under_approximation(const PartitionTreeSet<R>& pts, const FiniteGrid<R>& fg);

    
    template<class R> std::ostream& operator<<(std::ostream& os, const GridCell<R>& gc);
    template<class R> std::ostream& operator<<(std::ostream& os, const GridBlock<R>& gb);
    template<class R> std::ostream& operator<<(std::ostream& os, const GridCellListSet<R>& gcls);
    template<class R> std::ostream& operator<<(std::ostream& os, const GridMaskSet<R>& gms);
    
    
  }
}


#include "grid_set.inline.h"

#endif /* ARIADNE_GRID_SET_H */
