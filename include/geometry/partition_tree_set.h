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

#ifndef _ARIADNE_PARTITION_TREE_SET_H
#define _ARIADNE_PARTITION_TREE_SET_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/iterator.h"
#include "../base/tribool.h"

#include "../combinatoric/binary_word.h"
#include "../combinatoric/binary_tree.h"
#include "../combinatoric/subdivision_tree_set.h"

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
    template<class R, template<class> class BS> class ListSet;

          
    template<class R>
    class PartitionTreeIterator 
      : public boost::iterator_facade<PartitionTreeIterator<R>,
                                      PartitionTreeCell<R>,
                                      boost::forward_traversal_tag,
                                      PartitionTreeCell<R> >
    {
      friend class PartitionTree<R>;
     public:
      PartitionTreeIterator(const Rectangle<R>& bb, 
                            const Combinatoric::SubdivisionSequence& ss, 
                            Combinatoric::BinaryTree::const_iterator i)
        : _bounding_box(bb), _subdivisions(ss), _base(i) { }
     private:
      bool equal(const PartitionTreeIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Rectangle<R> _bounding_box;
      const Combinatoric::SubdivisionSequence _subdivisions;
      Combinatoric::BinaryTree::const_iterator _base;
    };

    template<class R>
    class PartitionTreeSetIterator 
      : public boost::iterator_facade<PartitionTreeSetIterator<R>,
                                      PartitionTreeCell<R>,
                                      boost::forward_traversal_tag,
                                      PartitionTreeCell<R> >
    {
      friend class PartitionTree<R>;
     public:
     private:
      bool equal(const PartitionTreeSetIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Rectangle<R> _bounding_box;
      const Combinatoric::SubdivisionSequence _subdivisions;
      Combinatoric::MaskedBinaryTree::const_iterator _base;
    };

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
      PartitionScheme(const Rectangle<R>& bb, const Combinatoric::SubdivisionSequence& sc)
        : _unit_box(bb), _subdivisions(sc) { }

      /*! \brief Equality. */
      bool operator==(const PartitionScheme<R>& pg) const {
        return _unit_box==pg._unit_box && _subdivisions==pg._subdivisions;
      }

      /*! \brief Inequality. */
      bool operator!=(const PartitionScheme<R>& pg) const {
        return !(*this==pg); 
      }
      
      /*! \brief The base unit box of the partition scheme. */
      const Rectangle<R>& unit_box() const { return _unit_box; }

      /*! \brief The sequence of subdivision coordinates. */
      const Combinatoric::SubdivisionSequence& subdivisions() const { return _subdivisions; }

      /*! \brief The underlying dimension of the partition scheme. */
      dimension_type dimension() const { return _subdivisions.dimension(); }
     private:
      Rectangle<R> _unit_box;
      Combinatoric::SubdivisionSequence _subdivisions;
    };



    /*!\ingroup BasicSet
     * \ingroup PartitionTree
     * \brief A rectangle defined on a partition tree.
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
      PartitionTreeCell(const Rectangle<R>& r, const Combinatoric::SubdivisionTreeCell& c)
        : _unit_box(r), _subdivision_cell(c)
      { assert(r.dimension()==c.dimension()); }

      /*!\brief Construct from a rectangle, the subdivision_coordinates and a binary word. */
      PartitionTreeCell(const Rectangle<R>& r, 
                        const Combinatoric::SubdivisionSequence& s, 
                        const Combinatoric::BinaryWord& w) 
        : _unit_box(r), _subdivision_cell(s,w) 
      { }

/*
      bool operator==(const PartitionTreeCell<R>& other) const {
        return this->_subdivision_cell==other._subdivision_cell 
          && this->_unit_box==other._unit_box;
      }
*/
      /*!\brief The cell in a unit box. */
      const Rectangle<R>& unit_box() const {
        return this->_unit_box; }

      /*!\brief The cell in a unit box. */
      const Combinatoric::SubdivisionTreeCell& subdivision_cell() const { 
        return this->_subdivision_cell; }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { 
        return this->_subdivision_cell.dimension(); }

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;
      
      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;
     private:
      const Rectangle<R> _unit_box;
      Combinatoric::SubdivisionTreeCell _subdivision_cell;
    };



    /*!\ingroup PartitionTree
     * \brief A tree structure following a PartitionScheme.
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
                             const Combinatoric::BinaryTree& t)
        : _unit_box(r), _subdivision_tree(s,t) { }

      /*! \brief Construct a tree based on a partition scheme and a binary tree. */
      explicit PartitionTree(const PartitionScheme<R>& ps, const Combinatoric::BinaryTree& t)
        : _unit_box(ps.unit_box()), _subdivision_tree(ps.subdivisions(),t) { }

      /*! \brief The unit box of the partition scheme. */
      const Rectangle<R>& unit_box() const { return _unit_box; }

      /*! \brief The underlying subdivision tree. */
      const Combinatoric::SubdivisionTree& subdivision_tree() const { return _subdivision_tree; }

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _subdivision_tree.dimension(); }

      /*! \brief The underlying bounding box. */
      const Combinatoric::SubdivisionSequence& subdivisions() const { return _subdivision_tree.subdivisions(); }

      /*! \brief The array describing the tree. */
      const Combinatoric::BinaryTree& binary_tree() const { return _subdivision_tree.binary_tree(); }

      /*! \brief The number of cells in the PartitionTree. */
      size_type size() const { return _subdivision_tree.size(); }

      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const { return  PartitionScheme<R>(unit_box(),subdivisions()); }
      
      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const { return const_iterator(_unit_box,_subdivision_tree.begin()); }
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const { return const_iterator(_unit_box,_subdivision_tree.end()); }
     private:
      Rectangle<R> _unit_box;
      Combinatoric::SubdivisionTree _subdivision_tree;
    };


    /*!\ingroup DenotableSet
     * \ingroup PartitionTree
     * \brief A denotable set on a partition grid, defined using a partition tree of cells.
     */
    template<class R>
    class PartitionTreeSet {
     public:
      //      typedef PartitionTreeSetIterator<R> iterator;
      //      typedef PartitionTreeSetIterator<R> const_iterator;
      typedef binary_constructor_iterator< Combinatoric::SubdivisionTreeSet::const_iterator,
                                           PartitionTreeCell<R>,
                                           Rectangle<R> > const_iterator;         
      typedef const_iterator iterator;

      /*! \brief Construct an empty set based on a partition scheme. */
      PartitionTreeSet(const PartitionScheme<R>& g)
        : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions()) 
      { }

      /*! \brief Construct an set based on a partition scheme, a binary tree and a mask. */
      PartitionTreeSet(const PartitionScheme<R>& g, const Combinatoric::BinaryTree& t, const BooleanArray& m)
        : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions(),t,m) 
      { }

      /*! \brief Construct a set based on a partition tree and a mask. */
      PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m)
        : _unit_box(t.unit_box()), _subdivision_set(t.subdivisions(),t.binary_tree(),m)
      { }

      /*! \brief Construct a set based on a bounding box, a subdivision sequence, a binary tree and a mask. */
      PartitionTreeSet(const Rectangle<R>& r, const Combinatoric::SubdivisionSequence& s, const Combinatoric::BinaryTree& t, const BooleanArray& m)
        : _unit_box(r), _subdivision_set(s,t,m)
      { }

      /*! \brief Convert from a GridMaskSet.
       *
       *  To ensure that the conversion is exact, and uses the same cells, the block covered by the mask must 
       *  have sides which are a power of two.
       */
      PartitionTreeSet(const GridMaskSet<R>& gms);

      /*! \brief Convert to a list of rectangles on a grid. */
      operator GridBlockListSet<R> () const;

      /*! \brief Convert to a ListSet of Rectangle. */
      operator ListSet<R,Rectangle> () const;

      /*! \brief A bounding box. */
      const Rectangle<R>& bounding_box() const { return _unit_box; }

      /*! \brief The unit box containing the partition tree set. */
      const Rectangle<R>& unit_box() const { return _unit_box; }

      /*! \brief The underlying subdivision set. */
      const Combinatoric::SubdivisionTreeSet& subdivision_set() const { return _subdivision_set; }

      /*! \brief The space dimension of the set. */
      size_type dimension() const { return _subdivision_set.dimension(); }

      /*! \brief The subdivision coordinates. */
      const Combinatoric::SubdivisionSequence& subdivisions() const { return _subdivision_set.subdivisions(); }

      /*! \brief The binary tree. */
      const Combinatoric::BinaryTree& binary_tree() const { return _subdivision_set.binary_tree(); }
      
      /*! \brief The mask. */
      const BooleanArray& mask() const { return _subdivision_set.mask(); }

      /*! \brief The number of cells in the partition tree. */
      size_type capacity() const { return _subdivision_set.capacity(); }

      /*! \brief The number of cells in the set. */
      size_type size() const { return _subdivision_set.size(); }

      /*! \brief The maximum depth in each coordinate. */
      SizeArray depths() const { return _subdivision_set.depths(); }

      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const { return _subdivision_set.depth(); }
      
      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const { 
        return PartitionScheme<R>(bounding_box(),subdivisions()); }
      
      /*! \brief The binary tree. */
      PartitionTree<R> partition_tree() const { 
        return PartitionTree<R>(bounding_box(),subdivisions(),binary_tree()); }
        
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { return const_iterator(_unit_box,_subdivision_set.begin()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_unit_box,_subdivision_set.end()); }
      
#ifdef DOXYGEN
      /*! \brief Compute an outer approximation to set \a s based on the partition scheme \a ps to depth a d. */
      friend template<class Set>
      PartitionTreeSet<R> outer_approximation(const Set& s, const PartitionScheme<R>& ps, const uint depth);

      /*! \brief Compute an inner approximation to set \a s based on the partition scheme \a ps to depth a d. */
      friend template<class Set>
      PartitionTreeSet<R> inner_approximation(const Set& s, const PartitionScheme<R>& ps, const uint depth);

      /*! \brief Compute an over approximation to set \a s based on the partition scheme \a ps to depth a d.
       *
       * This function may fail if the set \a s is not regular. 
       */
      friend template<class Set>
      PartitionTreeSet<R> over_approximation(const Set& s, const PartitionScheme<R>& ps, const uint depth);

      /*! \brief Compute an under approximation to set \a s based on the partition scheme \a ps to depth a d. */
      friend template<class Set>
      PartitionTreeSet<R> under_approximation(const Set& s, const PartitionScheme<R>& ps, const uint depth);
#endif
     private:
      Rectangle<R> _unit_box;
      Combinatoric::SubdivisionTreeSet _subdivision_set;
    };

    
    template<class R, class S>
    PartitionTreeSet<R> outer_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<class R, class S>
    PartitionTreeSet<R> inner_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<class R, class S>
    PartitionTreeSet<R> over_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<class R, class S>
    PartitionTreeSet<R> under_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
  }
}

#endif /* _ARIADNE_PARTITION_TREE_SET_H */
