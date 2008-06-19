/***************************************************************************
 *            grid_tree_set.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *
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

/*! \file grid_tree_set.h
 *  \brief Cuboidal grid trees.
 */

#ifndef ARIADNE_GRID_TREE_SET_H
#define ARIADNE_GRID_TREE_SET_H

#include <iosfwd>

#include "base/iterator.h"
#include "base/tribool.h"

#include "combinatoric/binary_word.h"
#include "combinatoric/binary_tree.h"

#include "geometry/declarations.h"
#include "geometry/set_interface.h"
#include "geometry/rectangle_expression.h"


namespace Ariadne {
  

    template<class R> class Grid;
    template<class R> class GridTreeCell;
    template<class R> class GridTreeSet;

    template<class R> class GridTreeSetIterator;

    /* External class declarations. */
    template<class R> class Box;
    template<class BS> class ListSet;


    /*!\ingroup BasicSet
     * \ingroup GridTree
     * \brief A rectangular cell in a grid tree.
     *
     * Defined as a SubdivisionTreeCell within a base cell given as a Box<R>.
     *
     * Satisfies the requirements of a RectangleExpression.
     */
    template<class R>
    class GridTreeCell
      : public RectangleExpression< GridTreeCell<R> >
    {
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real number used for the cell blounds. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by cells. */
      typedef Point<R> state_type;

      /*!\brief Construct from a rectangle, and a unit grid tree cell. */
      GridTreeCell(const Box<R>& r, const SubdivisionTreeCell& c);

      /*!\brief Construct from a rectangle, the subdivision_coordinates and a binary word. */
      GridTreeCell(const Box<R>& r, 
                        const SubdivisionSequence& s, 
                        const BinaryWord& w);

      /*!\brief The cell in a unit box. */
      const Box<R>& unit_box() const;

      /*!\brief The cell in a unit box. */
      const SubdivisionTreeCell& subdivision_cell() const;

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const;

      /*!\brief Tests if the set is empty. */
      tribool empty() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;
      
      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;

      /*! \brief An approximation to the volume of the set. */
      R volume() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      Box<R> _unit_box;
      SubdivisionTreeCell _subdivision_cell;
    };





    /*!\ingroup DenotableSet
     * \ingroup Grid
     * \brief A denotable set on a grid, defined using a tree of cells.
     *
     * Intersection (as an open set), union and set difference can be 
     * computed in time which is linear in the number elements of the grids
     * underlying the two sets.
     *
     * Testing inclusion of a cell in a %GridTreeSet, or adjoining a single
     * cell also takes unit time, which means that a LatticeMaskSet may be
     * preferable if these operations are common.
     *
     * Defined as a SubdivisionTree within a base cell given as a Box<R>.
     */
    template<class R>
    class GridTreeSet 
      : public SetInterface<R>
    {
     public:
      //      typedef GridTreeSetIterator<R> iterator;
      //      typedef GridTreeSetIterator<R> const_iterator;
      typedef binary_constructor_iterator< SubdivisionTreeSet::const_iterator,
                                           GridTreeCell<R>,
                                           Box<R> > const_iterator;         
      typedef const_iterator iterator;

      /*! \brief A tag describing the type of set. */
      typedef denotable_set_tag set_category;
      /*! \brief A tag describing the type of set. */
      typedef GridTreeCell<R> basic_set_type;


      /*! \brief Construct an empty set based on a grid. */
      GridTreeSet(const Grid<R>& g);

      /*! \brief Construct an set based on a grid scheme, a binary tree and a mask. */
      GridTreeSet(const Grid<R>& g, const BinaryTree& t, const BooleanArray& m);

      /*! \brief Convert from a GridMaskSet.
       *
       *  To ensure that the conversion is exact, and uses the same cells, the block covered by the mask must 
       *  have sides which are a power of two.
       */
      GridTreeSet(const GridMaskSet<R>& gms);

      //@{
      //! \name SetInterface methods
      /*! \brief Tests for disjointness with a Box. */
      virtual GridTreeSet<R>* clone() const;

      /*! \brief The space dimension of the set. */
      virtual dimension_type dimension() const;

      /*!\brief Checks if a denotable set includes a point. */
      virtual tribool contains(const Point<R>& p) const;

      /*! \brief Tests for superset of a Box. */ 
      virtual tribool superset(const Box<R>& r) const;

      /*! \brief Tests for intersection with a Box. */
      virtual tribool intersects(const Box<R>& r) const;

      /*! \brief Tests for disjointness with a Box. */
      virtual tribool disjoint(const Box<R>& r) const;

      /*! \brief Tests for subset of a Box. */
      virtual tribool subset(const Box<R>& r) const;

      /*!\brief A rectangle containing the grid cell. */
      virtual Box<R> bounding_box() const;
      //@}

      /*! \brief Convert to a ListSet of Box. */
      operator ListSet< Box<R> > () const;

      /*! \brief The underlying grid. */
      const Grid<R>& grid() const;

      /*! \brief The binary tree. */
      const BinaryTree& binary_tree() const;
      
      /*! \brief The mask. */
      const BooleanArray& mask() const;

      /*! \brief The number of cells in the tree. */
      size_type capacity() const;

      /*! \brief The number of cells in the set. */
      size_type size() const;

      /*! \brief The maximum depth in each coordinate. */
      SizeArray depths() const;

      /*! \brief The height of the top cell in the set. */
      size_type height() const;
      
      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const;
      
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const;
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const;
      
      /*!\brief Tests if the set is empty. */
      tribool empty() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*! \brief An approximation to the volume of the set. */
      R volume() const;
      
#ifdef DOXYGEN
      /*! \brief Compute an over approximation to a rectangle \a r based on the grid \a g. */
      friend template<class R>
      GridTreeCell<R> over_approximation(const Box<R>& r, const Grid<R>& g, const uint depth);

      /*! \brief Compute an outer approximation to set \a s based on the grid scheme \a ps to depth a d. */
      friend template<class SetInterface>
      GridTreeSet<R> outer_approximation(const SetInterface& s, const Grid<R>& ps, const uint depth);

      /*! \brief Compute an inner approximation to set \a s based on the grid scheme \a ps to depth a d. */
      friend template<class SetInterface>
      GridTreeSet<R> inner_approximation(const SetInterface& s, const Grid<R>& ps, const uint depth);
#endif

      /*! \brief Write a summary to an output stream. */
      std::ostream& summarize(std::ostream&) const;
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      static void _instantiate_geometry_operators();
     private:
      Grid<R> _grid;
      BinaryTree _tree;
      BooleanArray _mask;
    };

    
    template<class R>
    tribool contains(const GridTreeSet<R>&, const Point<R>&);

    template<class R>
    tribool disjoint(const GridTreeSet<R>&, const Box<R>&);

    template<class R>
    tribool subset(const GridTreeSet<R>&, const Box<R>&);

    template<class R>
    tribool subset(const Box<R>&, const GridTreeSet<R>&);

    
    template<class R>
    GridTreeCell<R> over_approximation(const Box<R>& r, const Grid<R>& g, const uint& d);
    
    template<class R, class S>
    GridTreeSet<R> outer_approximation(const S& s, const Grid<R>& g, const uint depth);
    
    template<class R, class S>
    GridTreeSet<R> inner_approximation(const S& s, const Grid<R>& ps, const uint depth);
    
    
    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const GridTreeCell<R>& ptc);
    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const GridTreeSet<R>& pts);
    
    
  }
}

#include "grid_tree_set.inline.h"

#endif /* ARIADNE_GRID_TREE_SET_H */
