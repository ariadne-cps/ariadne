/***************************************************************************
 *            partition_tree_set.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

/*! \file partition_tree_set.h
 *  \brief Cuboidal partition trees.
 */

#ifndef ARIADNE_PARTITION_TREE_SET_H
#define ARIADNE_PARTITION_TREE_SET_H

#include <iosfwd>

#include "../base/iterator.h"
#include "../base/tribool.h"

#include "../combinatoric/binary_word.h"
#include "../combinatoric/binary_tree.h"
#include "../combinatoric/subdivision_tree_set.h"

#include "../geometry/declarations.h"
#include "../geometry/set_interface.h"
#include "../geometry/rectangle_expression.h"


namespace Ariadne {
  namespace Geometry {

    template<class R> class PartitionScheme;
    template<class R> class PartitionTree;
    template<class R> class PartitionTreeCell;
    template<class R> class PartitionTreeSet;

    template<class R> class PartitionTreeSetIterator;

    template<class R> std::ostream& operator<<(std::ostream&, const PartitionScheme<R>&);
    template<class R> std::ostream& operator<<(std::ostream&, const PartitionTree<R>&);
    template<class R> std::ostream& operator<<(std::ostream&, const PartitionTreeCell<R>&);
    template<class R> std::ostream& operator<<(std::ostream&, const PartitionTreeSet<R>&);

    /* External class declarations. */
    template<class R> class Rectangle;
    template<class BS> class ListSet;

    template<class R> class PartitionTreeIterator;
    template<class R> class PartitionTreeSetIterator;

    /*!\ingroup PartitionTree
     * \brief A scheme for creating sets based on binary partition trees.
     */
    template<class R> 
    class PartitionScheme {
    public:
      /*! \brief The type of denotable real number used for the cell vertices. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by cells. */
      typedef Point<R> state_type;

      /*! \brief Construct from a bounding box \a bb with default subdivision coordinates. */
      PartitionScheme(const Rectangle<R>& bb);

      /*! \brief Construct from a bounding box \a bb and a sequence of subdivision coordinates \a sc. */
      PartitionScheme(const Rectangle<R>& bb, const Combinatoric::SubdivisionSequence& sc);

      /*! \brief Equality. */
      bool operator==(const PartitionScheme<R>& pg) const;

      /*! \brief Inequality. */
      bool operator!=(const PartitionScheme<R>& pg) const;
      
      /*! \brief The base unit box of the partition scheme. */
      const Rectangle<R>& unit_box() const;

      /*! \brief The sequence of subdivision coordinates. */
      const Combinatoric::SubdivisionSequence& subdivisions() const;

      /*! \brief The underlying dimension of the partition scheme. */
      dimension_type dimension() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      Rectangle<R> _unit_box;
      Combinatoric::SubdivisionSequence _subdivisions;
    };



    /*!\ingroup BasicSet
     * \ingroup PartitionTree
     * \brief A rectangular cell in a partition tree.
     *
     * Defined as a SubdivisionTreeCell within a base cell given as a Rectangle<R>.
     *
     * Satisfies the requirements of a RectangleExpression.
     */
    template<class R>
    class PartitionTreeCell
      : public RectangleExpression< PartitionTreeCell<R> >
    {
     public:
      /*! \brief The type of denotable real number used for the cell blounds. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by cells. */
      typedef Point<R> state_type;

      /*!\brief Construct from a rectangle, and a unit partition tree cell. */
      PartitionTreeCell(const Rectangle<R>& r, const Combinatoric::SubdivisionTreeCell& c);

      /*!\brief Construct from a rectangle, the subdivision_coordinates and a binary word. */
      PartitionTreeCell(const Rectangle<R>& r, 
                        const Combinatoric::SubdivisionSequence& s, 
                        const Combinatoric::BinaryWord& w);

      /*!\brief The cell in a unit box. */
      const Rectangle<R>& unit_box() const;

      /*!\brief The cell in a unit box. */
      const Combinatoric::SubdivisionTreeCell& subdivision_cell() const;

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
      const Rectangle<R> _unit_box;
      Combinatoric::SubdivisionTreeCell _subdivision_cell;
    };



    /*!\ingroup PartitionTree
     * \brief A tree structure following a PartitionScheme.
     *
     * Defined as a SubdivisionTree within a base cell given as a Rectangle<R>.
     */
    template<class R>
    class PartitionTree {
      friend class PartitionTreeSet<R>;
     public:
      typedef binary_constructor_iterator< Combinatoric::SubdivisionTree::const_iterator,
                                           PartitionTreeCell<R>,
                                           Rectangle<R> > const_iterator;
      typedef const_iterator iterator;

      /*! \brief The type of denotable real number used for the cell vertices. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by cells. */
      typedef Point<R> state_type;

      /*! \brief Construct a tree from a rectangle, a subdivision sequence and a binary tree. */
      explicit PartitionTree(const Rectangle<R>& r, 
                             const Combinatoric::SubdivisionSequence& s, 
                             const Combinatoric::BinaryTree& t);

      /*! \brief Construct a tree based on a partition scheme and a binary tree. */
      explicit PartitionTree(const PartitionScheme<R>& ps, const Combinatoric::BinaryTree& t);

      /*! \brief The unit box of the partition scheme. */
      const Rectangle<R>& unit_box() const;

      /*! \brief The underlying subdivision tree. */
      const Combinatoric::SubdivisionTree& subdivision_tree() const;

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const;

      /*! \brief The underlying bounding box. */
      const Combinatoric::SubdivisionSequence& subdivisions() const;

      /*! \brief The array describing the tree. */
      const Combinatoric::BinaryTree& binary_tree() const;

      /*! \brief The number of cells in the PartitionTree. */
      size_type size() const;

      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const;
      
      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const;
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      Rectangle<R> _unit_box;
      Combinatoric::SubdivisionTree _subdivision_tree;
    };


    /*!\ingroup DenotableSet
     * \ingroup PartitionTree
     * \brief A denotable set on a partition grid, defined using a partition tree of cells.
     *
     * Intersection (as an open set), union and set difference can be 
     * computed in time which is linear in the number elements of the partitions
     * underlying the two sets.
     *
     * Testing inclusion of a cell in a %PartitionTreeSet, or adjoining a single
     * cell also takes unit time, which means that a LatticeMaskSet may be
     * preferable if these operations are common.
     *
     * Defined as a SubdivisionTree within a base cell given as a Rectangle<R>.
     */
    template<class R>
    class PartitionTreeSet 
      : public SetInterface<R>
    {
     public:
      //      typedef PartitionTreeSetIterator<R> iterator;
      //      typedef PartitionTreeSetIterator<R> const_iterator;
      typedef binary_constructor_iterator< Combinatoric::SubdivisionTreeSet::const_iterator,
                                           PartitionTreeCell<R>,
                                           Rectangle<R> > const_iterator;         
      typedef const_iterator iterator;

      /*! \brief Construct an empty set based on a partition scheme. */
      PartitionTreeSet(const PartitionScheme<R>& g);

      /*! \brief Construct an set based on a partition scheme, a binary tree and a mask. */
      PartitionTreeSet(const PartitionScheme<R>& g, const Combinatoric::BinaryTree& t, const BooleanArray& m);

      /*! \brief Construct a set based on a partition tree and a mask. */
      PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m);

      /*! \brief Construct a set based on a bounding box, a subdivision sequence, a binary tree and a mask. */
      PartitionTreeSet(const Rectangle<R>& r, const Combinatoric::SubdivisionSequence& s, const Combinatoric::BinaryTree& t, const BooleanArray& m);

      /*! \brief Convert from a GridMaskSet.
       *
       *  To ensure that the conversion is exact, and uses the same cells, the block covered by the mask must 
       *  have sides which are a power of two.
       */
      PartitionTreeSet(const GridMaskSet<R>& gms);

      //@{
      //! \name SetInterface methods
      /*! \brief Tests for disjointness with a Rectangle. */
      virtual PartitionTreeSet<R>* clone() const;

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

      /*!\brief A rectangle containing the grid cell. */
      virtual Rectangle<R> bounding_box() const;
      //@}

      /*! \brief Convert to a ListSet of Rectangle. */
      operator ListSet< Rectangle<R> > () const;

      /*! \brief The unit box containing the partition tree set. */
      const Rectangle<R>& unit_box() const;

      /*! \brief The underlying subdivision set. */
      const Combinatoric::SubdivisionTreeSet& subdivision_set() const;

      /*! \brief The subdivision coordinates. */
      const Combinatoric::SubdivisionSequence& subdivisions() const;

      /*! \brief The binary tree. */
      const Combinatoric::BinaryTree& binary_tree() const;
      
      /*! \brief The mask. */
      const BooleanArray& mask() const;

      /*! \brief The number of cells in the partition tree. */
      size_type capacity() const;

      /*! \brief The number of cells in the set. */
      size_type size() const;

      /*! \brief The maximum depth in each coordinate. */
      SizeArray depths() const;

      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const;
      
      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const;
      
      /*! \brief The binary tree. */
      PartitionTree<R> partition_tree() const;
        
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
      /*! \brief Compute an outer approximation to set \a s based on the partition scheme \a ps to depth a d. */
      friend template<class SetInterface>
      PartitionTreeSet<R> outer_approximation(const SetInterface& s, const PartitionScheme<R>& ps, const uint depth);

      /*! \brief Compute an inner approximation to set \a s based on the partition scheme \a ps to depth a d. */
      friend template<class SetInterface>
      PartitionTreeSet<R> inner_approximation(const SetInterface& s, const PartitionScheme<R>& ps, const uint depth);

      /*! \brief Compute an over approximation to set \a s based on the partition scheme \a ps to depth a d.
       *
       * This function may fail if the set \a s is not regular. 
       */
      friend template<class SetInterface>
      PartitionTreeSet<R> over_approximation(const SetInterface& s, const PartitionScheme<R>& ps, const uint depth);

      /*! \brief Compute an under approximation to set \a s based on the partition scheme \a ps to depth a d. */
      friend template<class SetInterface>
      PartitionTreeSet<R> under_approximation(const SetInterface& s, const PartitionScheme<R>& ps, const uint depth);
#endif

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      static void _instantiate_geometry_operators();
     private:
      Rectangle<R> _unit_box;
      Combinatoric::SubdivisionTreeSet _subdivision_set;
    };

    
    template<class R>
    tribool contains(const PartitionTreeSet<R>&, const Point<R>&);

    template<class R>
    tribool disjoint(const PartitionTreeSet<R>&, const Rectangle<R>&);

    template<class R>
    tribool subset(const PartitionTreeSet<R>&, const Rectangle<R>&);

    template<class R>
    tribool subset(const Rectangle<R>&, const PartitionTreeSet<R>&);

    
    template<class R, class S>
    PartitionTreeSet<R> outer_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<class R, class S>
    PartitionTreeSet<R> inner_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<class R, class S>
    PartitionTreeSet<R> over_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<class R, class S>
    PartitionTreeSet<R> under_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    
    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const PartitionScheme<R>& ps);
    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const PartitionTreeCell<R>& ptc);
    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const PartitionTree<R>& pt);
    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const PartitionTreeSet<R>& pts);
    
    
  }
}

#include "partition_tree_set.inline.h"

#endif /* ARIADNE_PARTITION_TREE_SET_H */
