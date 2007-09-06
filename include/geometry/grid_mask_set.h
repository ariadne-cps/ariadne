/***************************************************************************
 *            grid_mask_set.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file grid_mask_set.h
 *  \brief A set on a grid defined by a mask.
 */

#ifndef ARIADNE_GRID_MASK_SET_H
#define ARIADNE_GRID_MASK_SET_H

#include <iosfwd>

#include <boost/iterator/iterator_adaptor.hpp>

#include "../base/array.h"
#include "../base/tribool.h"
#include "../base/pointer.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/declarations.h"
#include "../geometry/set_interface.h"
#include "../geometry/grid.h"
#include "../geometry/grid_cell.h"

/*TODO: Unify bounds in FiniteGrid, and make GridMaskSet use only finite grid bounds*/

namespace Ariadne {
  namespace Geometry {
      
    class denotable_set_tag;
    template<class Base, class Value> class GridSetIterator;

    
    
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
     *  To convert to a PartitionTreeSet, the number of cells used in each coordinate must be a power of two. 
     *  (This restriction should be lifted in a future version.)
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

      /*! \brief A tag describing the type of set. */
      typedef denotable_set_tag set_category;
      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      /*! \brief The type of basic set used in the list. */
      typedef GridCell<R> basic_set_type;
      /*! \brief The type of basic set contained by the denotable set. */
      typedef GridCell<R> value_type;
      /*! \brief The type of object returned by indexing. */

      /*!\brief Construct an empty set from a finite grid. (Deprecated)
       *
       * \deprecated Use GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b) instead.
       */
      explicit GridMaskSet(const FiniteGrid<R>& g);
     
      /*!\brief Construct a set from a finite grid and a mask. (Deprecated)
       *
       * \deprecated Use GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b, const BooleanArray& m) instead.
       */
      explicit GridMaskSet(const FiniteGrid<R>& g, const BooleanArray& m);

      /*!\brief Construct an empty set based on grid \a g such that the restricted grid covers bounding box \a bb as tightly as possible. */
      explicit GridMaskSet(const Grid<R>& g, const Rectangle<R>& bb);
     
      /*!\brief Construct an empty set based on grid \a g and with cells in the block \a b. */
      explicit GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b);
     
      /*!\brief Construct a set from a grid, a lattice rectangle and a mask. */
      explicit GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b, const BooleanArray& m);
     
      /*!\brief Construct a set from a grid, and a lattice mask set. */
      explicit GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeMaskSet& ms);
     
      /*!\brief Construct from a %GridCellListSet. Implicit conversion is disallowed since the grid bounds are synthesised automatically, which is usually not what is required. */
      explicit GridMaskSet(const GridCellListSet<R>& gcls);

      /*!\brief Copy constructor. */
      GridMaskSet(const GridMaskSet<R>& gms);

      /*!\brief Copy assignment. */
      GridMaskSet<R>& operator=(const GridMaskSet<R>& gms);

      //@{
      //! \name SetInterface methods
      /*! \brief Make a dynamically-allocated copy of the set. */
      virtual GridMaskSet<R>* clone() const;

      /*! \brief The space dimension of the set. */
      virtual dimension_type dimension() const;

      /*!\brief Checks if a denotable set includes a point. */
      virtual tribool contains(const Point<R>& p) const;

      /*! \brief Tests for superset of a Rectangle. */ 
      virtual tribool superset(const Rectangle<R>& r) const;

      /*! \brief Tests for intersection with a Rectangle. */
      virtual tribool intersects(const Rectangle<R>& r) const;

      /*! \brief Tests for disjointness with a Rectangle. */
      virtual tribool disjoint(const Rectangle<R>& r) const;

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

      /*! \brief The underlying lattice set. */
      const Combinatoric::LatticeMaskSet& lattice_set() const;

      /*! \brief The underlying finite_grid. */
      FiniteGrid<R> finite_grid() const;

      /*! \brief The bounds on elements in the set. */
      GridBlock<R> bounds() const;

      /*! \brief The extent of the elements in the grid. */
      Rectangle<R> extent() const;

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
      
      /*! \brief Returns true if the set is bounded. */
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

      /*! \brief Makes the set unbounded. */
      void adjoin_unbounded_cell();
    
      /*! \brief Removes a cell from the set. */
      void remove(const GridCell<R>& gc);

      /*! \brief Removes a list of cells from the set. */
      void remove(const GridCellListSet<R>& gcls);

      /*! \brief Removes a GridMaskSet from the set. */
      void remove(const GridMaskSet<R>& gms);

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

      /*! \brief Adjoins an over-approximation of the rectangle \a r. */
      void adjoin_over_approximation(const Rectangle<R>& r);

      /*! \brief Adjoins an under-approximation of the rectangle \a r. */
      void adjoin_under_approximation(const Rectangle<R>& r);

      /*! \brief Adjoins an outer-approximation of the basic set \a bs. */
      template<class BS>
      void adjoin_outer_approximation(const BS& bs);

      /*! \brief Adjoins an inner-approximation of the basic set \a bs. */
      template<class BS>
      void adjoin_inner_approximation(const BS& bs);

      /*! \brief Adjoins an outer-approximation of the set \a s. */
      void adjoin_outer_approximation(const SetInterface<R>& s);

      /*! \brief Adjoins an inner-approximation of the set \a s. */
      void adjoin_inner_approximation(const SetInterface<R>& s);

      /*! \brief Restricts to an outer-approximation of the set \a s. */
      void restrict_outer_approximation(const SetInterface<R>& s);

      /*! \brief Restricts to an inner-approximation of the set \a s. */
      void restrict_inner_approximation(const SetInterface<R>& s);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;

#ifdef DOXYGEN
      friend tribool overlap<> (const GridBlock<R>&, const GridMaskSet<R>&);
      friend tribool overlap<> (const GridMaskSet<R>&, const GridBlock<R>&);
      friend tribool overlap<> (const GridCellListSet<R>&, const GridMaskSet<R>&);
      friend tribool overlap<> (const GridMaskSet<R>&, const GridCellListSet<R>&);
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
      Base::shared_ptr< Grid<R> > _grid_ptr;
      Combinatoric::LatticeMaskSet _lattice_set;
    };
    
    
    template<class R> tribool disjoint(const Rectangle<R>&, const GridMaskSet<R>&);
    template<class R> tribool disjoint(const GridMaskSet<R>& gm, const Rectangle<R>& r);
    template<class R> tribool subset(const Rectangle<R>&, const GridMaskSet<R>&);
    template<class R> tribool subset(const GridMaskSet<R>&, const Rectangle<R>&);
   
    template<class R> tribool overlap(const GridBlock<R>&, const GridMaskSet<R>&);
    template<class R> tribool overlap(const GridMaskSet<R>&, const GridBlock<R>&);
    template<class R> tribool overlap(const GridCellListSet<R>&, const GridMaskSet<R>&);
    template<class R> tribool overlap(const GridMaskSet<R>&, const GridCellListSet<R>&);
    template<class R> tribool overlap(const GridMaskSet<R>&, const GridMaskSet<R>&);
    
    template<class R> tribool subset(const GridCell<R>&, const GridMaskSet<R>&);
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

   
 
    template<class R, class BS> GridMaskSet<R> outer_approximation(const ListSet<BS>& ls, const FiniteGrid<R>& fg);
 
    template<class R> GridMaskSet<R> outer_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg);
    template<class R> GridMaskSet<R> inner_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg);

    template<class R> ListSet< Rectangle<R> > lower_approximation(const SetInterface<R>& set, const Grid<R>& fg);
    template<class R> ListSet< Rectangle<R> > lower_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg);
 

    template<class R> std::ostream& operator<<(std::ostream& os, const GridMaskSet<R>& gms);
    
    
  }
}


#include "grid_mask_set.inline.h"

#endif /* ARIADNE_GRID_MASK_SET_H */
