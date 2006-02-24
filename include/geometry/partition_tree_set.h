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
#include <algorithm>
#include <vector>
#include <iostream>

#include "base/basic_type.h"
#include "base/binary_word.h"
#include "base/binary_tree.h"
#include "base/utility.h"

#include "geometry/geometry_declarations.h"
#include "geometry/unit_partition_tree_set.h"

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

          
    /*!\brief A partition grid of rectangles in Euclidean space.
     */
    template<typename R> class PartitionScheme {
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

      /*! \brief The underlying dimension of the partition scheme. */
      dimension_type dimension() const { return _bounding_box.dimension(); }

      /*! \brief The outer bounding box of the grid. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }
      /*! \brief The sequence of subdivision coordinates. */
      const SubdivisionSequence& subdivisions() const { return _subdivisions; }
     private:
      Rectangle<R> _bounding_box;
      SubdivisionSequence _subdivisions;
    };

    template<typename R>
    class PartitionTreeIterator {
      typedef PartitionTreeIterator<R> Self;
      friend class PartitionTree<R>;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef PartitionTreeCell<R> value_type;
      typedef PartitionTreeCell<R> reference;
      typedef const PartitionTreeCell<R>* pointer;
      typedef int difference_type;
     private:
      PartitionTreeIterator(const Rectangle<R>& bb, const SubdivisionSequence& ss, BinaryTree::const_iterator i)
        : _bounding_box(bb), _subdivisions(ss), _iter(i) { }
     public:
      bool operator==(const PartitionTreeIterator<R>& other) const { 
        return this->equal(other); }
      bool operator!=(const PartitionTreeIterator<R>& other) const { 
        return !this->equal(other); }
      PartitionTreeIterator<R>& operator++() { 
        this->increment(); return *this; }
      PartitionTreeIterator<R> operator++(int) { 
        PartitionTreeIterator<R> tmp(*this); this->increment(); return tmp; }
      PartitionTreeCell<R> operator*() const { 
        return this->dereference(); }
     private:
      bool equal(const PartitionTreeIterator<R>& other) const {
        return this->_iter == other._iter; }
      void increment() { ++_iter; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_iter); }
     private:
      const Rectangle<R> _bounding_box;
      const SubdivisionSequence _subdivisions;
      BinaryTreeIterator _iter;
    };


    /*! \brief A tree structure following a PartitionScheme.
     */
    template<typename R>
    class PartitionTree {
      friend class PartitionTreeSet<R>;
     public:
      typedef PartitionTreeIterator<R> iterator;
      typedef PartitionTreeIterator<R> const_iterator;

      /*! \brief Construct a tree from a rectangle, a subdivision sequence and a binary tree. */
      explicit PartitionTree(const Rectangle<R>& r, const SubdivisionSequence& s, const BinaryTree& t)
        : _bounding_box(r), _subdivisions(s), _tree(t) { }

      /*! \brief Construct a tree based on a partition scheme and a binary tree. */
      explicit PartitionTree(const PartitionScheme<R>& g, const BinaryTree& t)
        : _bounding_box(g.bounding_box()), _subdivisions(g.subdivisions()),_tree(t) { }

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _bounding_box.dimension(); }

      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const { return  PartitionScheme<R>(_bounding_box,_subdivisions); }
      
      /*! \brief The underlying bounding box. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }

      /*! \brief The underlying bounding box. */
      const SubdivisionSequence& subdivisions() const { return _subdivisions; }

      /*! \brief The array describing the tree. */
      const BinaryTree& binary_tree() const { return _tree; }

      /*! \brief The number of cells in the PartitionTree. */
      size_type size() const { return _tree.size(); }

      /*! \brief The number of cells in the PartitionTree. */
      size_type capacity() const { return _tree.size(); }

      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const { return const_iterator(_bounding_box,_subdivisions,_tree.begin()); }
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const { return const_iterator(_bounding_box,_subdivisions,_tree.end()); }
     private:
      Rectangle<R> _bounding_box;
      SubdivisionSequence _subdivisions;
      BinaryTree _tree;
    };



    /*! \brief A rectangle defined on a partition tree.
     */
    template<typename R>
    class PartitionTreeCell {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
      typedef UnitPartitionTreeCell::dyadic_type dyadic_type;

      /*!\brief Construct from a rectangle, and a unit partition tree cell. */
      PartitionTreeCell(const Rectangle<R>& r, const UnitPartitionTreeCell& c)
        : _bounding_box(r), _cell(c) { assert(r.dimension()==c.dimension()); }

      /*!\brief Construct from a rectangle, the subdivision_coordinates and a binary word. */
      PartitionTreeCell(const Rectangle<R>& r, const SubdivisionSequence& s, 
                        const BinaryWord& w) 
        : _bounding_box(r), _cell(r.dimension(),s,w) { }

      bool operator==(const PartitionTreeCell<R>& other) const {
        return this->_bounding_box==other._bounding_box && this->_cell==other._cell;
      }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return this->_bounding_box.dimension(); }

      /*!\brief The bounding box. */
      const Rectangle<R>& bounding_box() const { return this->_bounding_box; }

      /*!\brief The sequence of subdivisions used to create the cell. */
      const SubdivisionSequence& subdivisions() const { return this->_cell.subdivisions(); }
      
      /*!\brief The binary word defining the cell. */
      BinaryWord word() const { return this->_cell.word(); }
      
      /*!\brief The subdivision coordinates. */
      const UnitPartitionTreeCell& unit_cell() const { return _cell; }

      /*!\brief The bounding box. */
      Rectangle<dyadic_type> bounds() const { return this->_cell.bounds(); }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      const Rectangle<R> _bounding_box;
      UnitPartitionTreeCell _cell;
    };


    template<typename R>
    class PartitionTreeSetIterator {
      typedef PartitionTreeSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef PartitionTreeCell<R> value_type;
      typedef PartitionTreeCell<R> reference;
      typedef const PartitionTreeCell<R>* pointer;
      typedef int difference_type;
     public:
      PartitionTreeSetIterator(const Rectangle<R>& bb, const UnitPartitionTreeSetIterator& i)
        : _bounding_box(bb), _iter(i) { }
      PartitionTreeSetIterator(const PartitionTreeSetIterator<R>& other)
        : _bounding_box(other._bounding_box), _iter(other._iter) { }
        
      bool operator==(const PartitionTreeSetIterator<R>& other) { return this->equal(other); }
      bool operator!=(const PartitionTreeSetIterator<R>& other) { return !this->equal(other); }
      PartitionTreeCell<R> operator*() const { return this->dereference(); }
      Self& operator++() { this->increment(); return *this; }
      Self operator++(int) { PartitionTreeSetIterator<R> tmp(*this); this->increment(); return tmp; }
     private:
      bool equal(const PartitionTreeSetIterator<R>& other) {
        return this->_iter==other._iter && this->_bounding_box==other._bounding_box; }
      void increment() { ++_iter; }
      PartitionTreeCell<R> dereference() const { return PartitionTreeCell<R>(_bounding_box,*_iter); }
     private:
      const Rectangle<R>& _bounding_box;
      UnitPartitionTreeSetIterator _iter;
     };


    /*! \brief A denotable set on a partition grid, defined using a partition tree of cells.
     */
    template<typename R>
    class PartitionTreeSet {
     public:
      typedef PartitionTreeSetIterator<R> iterator;
      typedef PartitionTreeSetIterator<R> const_iterator;

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

      /*! \brief Convert from a GridMaskSet. */
      PartitionTreeSet(const GridMaskSet<R>& gms);

      /*! \brief Convert to a list of rectangles on a grid. */
      //operator GridRectangleListSet<R> () const;

      /*! \brief Convert to a ListSet of Rectangle. */
      operator ListSet<R,Rectangle> () const;

      /*! \brief The space dimension of the set. */
      size_type dimension() const { return _bounding_box.dimension(); }

      /*! \brief The underlying PartitionScheme. */
      PartitionScheme<R> scheme() const { 
        return PartitionScheme<R>(_bounding_box,_unit_set.subdivisions()); }
      
      /*! \brief The binary tree. */
      PartitionTree<R> partition_tree() const { 
        return PartitionTree<R>(_bounding_box,_unit_set.subdivisions(),_unit_set.tree()); }
        
      /*! \brief The bounding box. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }

      /*! \brief The subdivision coordinates. */
      const SubdivisionSequence& subdivisions() const { return _unit_set.subdivisions(); }

      /*! \brief The binary tree. */
      const BinaryTree& binary_tree() const { return _unit_set.tree(); }
      
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
      
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { return const_iterator(_bounding_box,_unit_set.begin()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_bounding_box,_unit_set.end()); }
     private:
      Rectangle<R> _bounding_box;
      UnitPartitionTreeSet _unit_set;
    };

  }
}

#endif /* _ARIADNE_PARTITION_TREE_SET_H */
