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

#include "../base/binary_word.h"
#include "../base/binary_tree.h"
#include "../base/iterator.h"

#include "../geometry/subdivision_tree_set.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R> class PartitionScheme;
    template<typename R> class PartitionTree;
    template<typename R> class PartitionTreeCell;
    template<typename R> class PartitionTreeSet;

    template<typename R> class PartitionTreeSetIterator;

    template<typename R> std::ostream& operator<<(std::ostream&, const PartitionScheme<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const PartitionTree<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const PartitionTreeCell<R>&);
    template<typename R> std::ostream& operator<<(std::ostream&, const PartitionTreeSet<R>&);

    /* External class declarations. */
    template<typename R> class Rectangle;
    template<typename R, template<typename> class BS> class ListSet;

          
    template<typename R>
    class PartitionTreeIterator 
      : public boost::iterator_facade<PartitionTreeIterator<R>,
                                      PartitionTreeCell<R>,
                                      boost::forward_traversal_tag,
                                      PartitionTreeCell<R> >
    {
      friend class PartitionTree<R>;
     public:
      PartitionTreeIterator(const Rectangle<R>& bb, const SubdivisionSequence& ss, BinaryTree::const_iterator i)
        : _bounding_box(bb), _subdivisions(ss), _base(i) { }
     private:
      bool equal(const PartitionTreeIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Rectangle<R> _bounding_box;
      const SubdivisionSequence _subdivisions;
      BinaryTree::const_iterator _base;
    };

    template<typename R>
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
      const SubdivisionSequence _subdivisions;
      MaskedBinaryTree::const_iterator _base;
    };

    /*! \brief A partition grid of rectangles in Euclidean space.
     *  \ingroup PartitionTree
     */
    template<typename R> 
    class PartitionScheme {
      typedef R real_type;
    public:
      /*! \brief Construct from a bounding box \a bb with default subdivision coordinates. */
      PartitionScheme(const Rectangle<R>& bb);

      /*! \brief Construct from a bounding box \a bb and a sequence of subdivision coordinates \a sc. */
      PartitionScheme(const Rectangle<R>& bb, const SubdivisionSequence& sc)
        : _bounding_box(bb), _subdivisions(sc) { }

      /*! \brief Equality. */
      bool operator==(const PartitionScheme<R>& pg) const {
        return _bounding_box==pg._bounding_box && _subdivisions==pg._subdivisions;
      }

      /*! \brief The outer bounding box of the grid. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }

      /*! \brief The sequence of subdivision coordinates. */
      const SubdivisionSequence& subdivisions() const { return _subdivisions; }

      /*! \brief The underlying dimension of the partition scheme. */
      dimension_type dimension() const { return _subdivisions.dimension(); }
     private:
      Rectangle<R> _bounding_box;
      SubdivisionSequence _subdivisions;
    };



    /*! \brief A rectangle defined on a partition tree.
     *  \ingroup BasicSet
     *  \ingroup PartitionTree
     */
    template<typename R>
    class PartitionTreeCell {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
      typedef SubdivisionTreeCell::dyadic_type dyadic_type;

      /*!\brief Construct from a rectangle, and a unit partition tree cell. */
      PartitionTreeCell(const Rectangle<R>& r, const SubdivisionTreeCell& c)
        : _bounding_box(r), _unit_cell(c)
      { assert(r.dimension()==c.dimension()); }

      /*!\brief Construct from a rectangle, the subdivision_coordinates and a binary word. */
      PartitionTreeCell(const Rectangle<R>& r, 
                        const SubdivisionSequence& s, 
                        const BinaryWord& w) 
        : _bounding_box(r), _unit_cell(s,w) 
      { }

      bool operator==(const PartitionTreeCell<R>& other) const {
        return this->_unit_cell==other._unit_cell 
          && this->_bounding_box==other._bounding_box;
      }

      /*!\brief The cell in a unit box. */
      const Rectangle<R>& bounding_box() const {
        return this->_bounding_box; }

      /*!\brief The cell in a unit box. */
      const SubdivisionTreeCell& unit_cell() const { 
        return this->_unit_cell; }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { 
        return this->_unit_cell.dimension(); }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Rectangle<R> _bounding_box;
      SubdivisionTreeCell _unit_cell;
    };



    /*! \brief A tree structure following a PartitionScheme.
     *  \ingroup PartitionTree
     */
    template<typename R>
    class PartitionTree {
      friend class PartitionTreeSet<R>;
     public:
      //      typedef PartitionTreeIterator<R> iterator;
      //      typedef PartitionTreeIterator<R> const_iterator;
      typedef binary_constructor_iterator< SubdivisionTree::const_iterator,
                                           PartitionTreeCell<R>,
                                           Rectangle<R> > const_iterator;
      typedef const_iterator iterator;

      /*! \brief Construct a tree from a rectangle, a subdivision sequence and a binary tree. */
      explicit PartitionTree(const Rectangle<R>& r, const SubdivisionSequence& s, const BinaryTree& t)
        : _bounding_box(r), _unit_tree(s,t) { }

      /*! \brief Construct a tree based on a partition scheme and a binary tree. */
      explicit PartitionTree(const PartitionScheme<R>& ps, const BinaryTree& t)
        : _bounding_box(ps.bounding_box()), _unit_tree(ps.subdivisions(),t) { }

      /*! \brief The underlying bounding box. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }

      /*! \brief The underlying bounding box. */
      const SubdivisionTree& unit_tree() const { return _unit_tree; }

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _unit_tree.dimension(); }

      /*! \brief The underlying bounding box. */
      const SubdivisionSequence& subdivisions() const { return _unit_tree.subdivisions(); }

      /*! \brief The array describing the tree. */
      const BinaryTree& binary_tree() const { return _unit_tree.binary_tree(); }

      /*! \brief The number of cells in the PartitionTree. */
      size_type size() const { return _unit_tree.size(); }

      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const { return  PartitionScheme<R>(bounding_box(),subdivisions()); }
      
      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const { return const_iterator(_bounding_box,_unit_tree.begin()); }
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const { return const_iterator(_bounding_box,_unit_tree.end()); }
     private:
      Rectangle<R> _bounding_box;
      SubdivisionTree _unit_tree;
    };


    /*! \brief A denotable set on a partition grid, defined using a partition tree of cells.
     *  \ingroup DenotableSet
     *  \ingroup PartitionTree
     */
    template<typename R>
    class PartitionTreeSet {
     public:
      //      typedef PartitionTreeSetIterator<R> iterator;
      //      typedef PartitionTreeSetIterator<R> const_iterator;
      typedef binary_constructor_iterator< SubdivisionTreeSet::const_iterator,
                                           PartitionTreeCell<R>,
                                           Rectangle<R> > const_iterator;         
      typedef const_iterator iterator;

      /*! \brief Construct an empty set based on a PartitionScheme. */
      PartitionTreeSet(const PartitionScheme<R>& g)
        : _bounding_box(g.bounding_box()), _unit_set(g.subdivisions()) 
      { }

      /*! \brief Construct an set based on a PartitionScheme, a binary tree and a mask. */
      PartitionTreeSet(const PartitionScheme<R>& g, const BinaryTree& t, const BooleanArray& m)
        : _bounding_box(g.bounding_box()), _unit_set(g.subdivisions(),t,m) 
      { }

      /*! \brief Construct a set based on a partition tree and a mask. */
      PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m)
        : _bounding_box(t.bounding_box()), _unit_set(t.subdivisions(),t.binary_tree(),m)
      { }

      /*! \brief Construct a set based on a bounding box, a subdivision sequence, a binary tree and a mask. */
      PartitionTreeSet(const Rectangle<R>& r, const SubdivisionSequence& s, const BinaryTree& t, const BooleanArray& m)
        : _bounding_box(r), _unit_set(s,t,m)
      { }

      /*! \brief Convert from a GridMaskSet.
       *
       *  To ensure that the conversion is exact, and uses the same cells, the lattice rectangle must 
       *  have sides which are a power of two.
       */
      PartitionTreeSet(const GridMaskSet<R>& gms);

      /*! \brief Convert to a list of rectangles on a grid. */
      //operator GridRectangleListSet<R> () const;

      /*! \brief Convert to a ListSet of Rectangle. */
      operator ListSet<R,Rectangle> () const;

      /*! \brief The bounding box. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }

      /*! \brief The bounding box. */
      const SubdivisionTreeSet& unit_set() const { return _unit_set; }

      /*! \brief The space dimension of the set. */
      size_type dimension() const { return _unit_set.dimension(); }

      /*! \brief The subdivision coordinates. */
      const SubdivisionSequence& subdivisions() const { return _unit_set.subdivisions(); }

      /*! \brief The binary tree. */
      const BinaryTree& binary_tree() const { return _unit_set.binary_tree(); }
      
      /*! \brief The mask. */
      const BooleanArray& mask() const { return _unit_set.mask(); }

      /*! \brief The number of cells in the partition tree. */
      size_type capacity() const { return _unit_set.capacity(); }

      /*! \brief The number of cells in the set. */
      size_type size() const { return _unit_set.size(); }

      /*! \brief The maximum depth in each coordinate. */
      SizeArray depths() const { return _unit_set.depths(); }

      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const { return _unit_set.depth(); }
      
      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const { 
        return PartitionScheme<R>(bounding_box(),subdivisions()); }
      
      /*! \brief The binary tree. */
      PartitionTree<R> partition_tree() const { 
        return PartitionTree<R>(bounding_box(),subdivisions(),binary_tree()); }
        
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { return const_iterator(_bounding_box,_unit_set.begin()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_bounding_box,_unit_set.end()); }
      
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
      Rectangle<R> _bounding_box;
      SubdivisionTreeSet _unit_set;
    };

    
    template<typename R, class S>
    PartitionTreeSet<R> outer_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<typename R, class S>
    PartitionTreeSet<R> inner_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<typename R, class S>
    PartitionTreeSet<R> over_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
    template<typename R, class S>
    PartitionTreeSet<R> under_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth);
    
  }
}

#endif /* _ARIADNE_PARTITION_TREE_SET_H */
