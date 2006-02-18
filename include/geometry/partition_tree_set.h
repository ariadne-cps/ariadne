/***************************************************************************
 *            partition_tree_set.cc
 *
 *  1 July 2006
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

    SubdivisionSequence default_subdivision_coordinates(dimension_type n);

    /* External class declarations. */
    template<typename R> class Rectangle;
    template<typename R, template<typename> class BS> class ListSet;

          
    /*!\brief A partition grid of rectangles in Euclidean space.
     */
    template<typename R> class PartitionScheme {
      typedef R real_type;
    public:
      /*! \brief Construct from a bounding box \a bb with default subdivision coordinates. */
      PartitionScheme(const Rectangle<R>& bb)
        : _bounding_box(bb), 
          _subdivision_coordinates(default_subdivision_coordinates(bb.dimension())) 
      { }

      /*! \brief Construct from a bounding box \a bb and a sequence of subdivision coordinates \a sc. */
      PartitionScheme(const Rectangle<R>& bb, const SubdivisionSequence& sc)
        : _bounding_box(bb), _subdivision_coordinates(sc) { }

      /*! \brief Equality. */
      bool operator==(const PartitionScheme<R>& pg) const {
        return _bounding_box==pg._bounding_box && _subdivision_coordinates==pg._subdivision_coordinates;
      }

      /*! \brief The underlying dimension of the partition scheme. */
      dimension_type dimension() const { return _bounding_box.dimension(); }

      /*! \brief The outer bounding box of the grid. */
      const Rectangle<R>& bounding_box() const { return _bounding_box; }
      /*! \brief The sequence of subdivision coordinates. */
      const SubdivisionSequence& subdivision_coordinates() const { return _subdivision_coordinates; }

      /*! \brief The coordinate in which the \a n th subdivision is performed. */
      dimension_type subdivision_coordinate(size_type n) const { return _subdivision_coordinates[n]; }
     private:
      Rectangle<R> _bounding_box;
      SubdivisionSequence _subdivision_coordinates;
    };



    /*! \brief A tree structure following a PartitionScheme.
     */
    template<typename R>
    class PartitionTree {
      friend class PartitionTreeSet<R>;
     public:
      class const_iterator {
        typedef const_iterator Self;
        friend class PartitionTree<R>;
       private:
        const_iterator(const PartitionScheme<R>& g, BinaryTree::const_iterator i) : _scheme(g), _iter(i) { }
       public:
        bool operator==(const Self& other) const { return this->_iter == other._iter; }
        bool operator!=(const Self& other) const { return this->_iter != other._iter; }
        const_iterator& operator++() { ++_iter; return *this; }
        PartitionTreeCell<R> operator*() const { return PartitionTreeCell<R>(_scheme,*_iter); }
       private:
        const PartitionScheme<R>& _scheme;
        BinaryTreeIterator _iter;
      };

      /*! \brief Construct a tree from a partition scheme and a binary tree. */
      explicit PartitionTree(const PartitionScheme<R>& g, const BinaryTree& t)
        : _scheme(g), _tree(t) { }

      /*! \brief Construct a tree from a rectangle, a subdivision sequence and a binary tree. */
      explicit PartitionTree(const Rectangle<R>& r, const SubdivisionSequence& s, const BinaryTree& t)
        : _scheme(r,s), _tree(t) { }

      /*! \brief Construct a tree from a rectangle and a binary tree. */
      explicit PartitionTree(const Rectangle<R>& r, const BinaryTree& t)
        : _scheme(r), _tree(t) { }

      /*! \brief Construct a one-cell tree. */
      explicit PartitionTree(const PartitionScheme<R>& g)
        : _scheme(g), _tree() { }

      /*! \brief Construct a one-cell tree. */
      explicit PartitionTree(const Rectangle<R>& bb, const SubdivisionSequence& ss)
        : _scheme(PartitionScheme<R>(bb,ss)), _tree() { }

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _scheme.dimension(); }

      /*! \brief The underlying PartitionScheme. */
      const PartitionScheme<R>& scheme() const { return _scheme; }

      /*! \brief The underlying bounding box. */
      const Rectangle<R>& bounding_box() const { return _scheme.bounding_box(); }

      /*! \brief The underlying bounding box. */
      const SubdivisionSequence& subdivision_coordinates() const { return _scheme.subdivision_coordinates(); }

      /*! \brief The coordinate in which the \a n th subdivision is performed. */
      dimension_type subdivision_coordinate(size_type n) const { return _scheme.subdivision_coordinate(n); }

      /*! \brief The array describing the tree. */
      const BinaryTree& tree() const { return _tree; }

      /*! \brief The number of cells in the PartitionTree. */
      size_type capacity() const { return _tree.size(); }

      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const { return const_iterator(_scheme,_tree.begin()); }
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const { return const_iterator(_scheme,_tree.end()); }
     private:
      PartitionScheme<R> _scheme;
      BinaryTree _tree;
    };



    /*! \brief A rectangle defined on a partition tree.
     */
    template<typename R>
    class PartitionTreeCell {
     public:
      typedef R real_type;
      typedef Point<R> state_type;

      /*!\brief Construct from a partition scheme and a binary word. */
      PartitionTreeCell(const PartitionScheme<R>& p, const BinaryWord& w)
        : _scheme(p), _word(w) { }

      /*!\brief Construct from a rectangle, the subdivision_coordinates and a binary word. */
      PartitionTreeCell(const Rectangle<R>& r, const SubdivisionSequence& s, 
                        const BinaryWord& w) 
        : _scheme(r,s), _word(w) { }

      bool operator==(const PartitionTreeCell<R>& ptc) const {
        return _scheme==ptc._scheme && _word==ptc._word;
      }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return _scheme.dimension(); }

      /*!\brief The underlying partition scheme. */
      const PartitionScheme<R>& scheme() const { return _scheme; }

      /*!\brief The bounding box. */
      const Rectangle<R>& bounding_box() const { return _scheme.bounding_box(); }

      /*!\brief The subdivision coordinates. */
      const SubdivisionSequence& subdivision_coordinates() const { return _scheme.subdivision_coordinates(); }

      /*!\brief The subdivision coordinates. */
      const BinaryWord& word() const { return _word; }

      /*!\brief The dimension of the cell. */
      BinaryWord index() const { return _word; }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      PartitionScheme<R> _scheme;
      BinaryWord _word;
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
      PartitionTreeSetIterator(const PartitionScheme<R>& ps, const BinarySubtreeIterator& i)
        : _bounding_box(ps.bounding_box()), _subdivision_coordinates(ps.subdivision_coordinates()), _iter(i) { }
      PartitionTreeSetIterator(const PartitionTreeSetIterator<R>& ptsi)
        : _bounding_box(ptsi._bounding_box), _subdivision_coordinates(ptsi._subdivision_coordinates), _iter(ptsi._iter) { }
      //FIXME: Check equality of sets
      bool operator==(const PartitionTreeSetIterator<R>& other) {
        return (this->_iter==other._iter); }
      bool operator!=(const PartitionTreeSetIterator<R>& other) { return !(*this==other); }
      PartitionTreeCell<R> operator*() const { return PartitionTreeCell<R>(_bounding_box,_subdivision_coordinates,*_iter); }
      Self& operator++() { ++_iter; return *this; }
      Self operator++(int) { PartitionTreeSetIterator<R> tmp(*this); ++(*this); return tmp; }
     private:
      const Rectangle<R>& _bounding_box;
      const SubdivisionSequence& _subdivision_coordinates;
      BinarySubtreeIterator _iter;
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
        : _ptree(g), _mask(1) { _mask[0]=false; }

      /*! \brief Construct a set based on a PartitionScheme, a binary tree and a mask. */
      PartitionTreeSet(const PartitionScheme<R>& g, const BinaryTree& t, const BooleanArray& m)
        : _ptree(g,t), _mask(m)
      {
        assert(_ptree.capacity() == _mask.size());
      }

      /*! \brief Construct a set based on a PartitionTree and a mask. */
      PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m)
        : _ptree(t), _mask(m)
      {
        assert(_ptree.capacity() == _mask.size());
      }

      /*! \brief Construct a set based on a bounding box, a subdivision sequence, a binary tree and a mask. */
      PartitionTreeSet(const Rectangle<R>& r, const SubdivisionSequence& s, const BinaryTree& t, const BooleanArray& m)
        : _ptree(r,s,t), _mask(m)
      {
        assert(_ptree.capacity() == _mask.size());
      }

      /*! \brief Convert from a GridMaskSet. */
      PartitionTreeSet(const GridMaskSet<R>& gms);

      /*! \brief Convert to a list of rectangles on a grid. */
      //operator GridRectangleListSet<R> () const;

      /*! \brief Convert to a ListSet of Rectangle. */
      operator ListSet<R,Rectangle> () const;

      /*! \brief The space dimension of the set. */
      size_type dimension() const { return scheme().dimension(); }

      /*! \brief The underlying PartitionScheme. */
      const PartitionScheme<R>& scheme() const { return _ptree.scheme(); }

      /*! \brief The bounding box. */
      const Rectangle<R>& bounding_box() const { return _ptree.bounding_box(); }

      /*! \brief The subdivision coordinates. */
      const SubdivisionSequence& subdivision_coordinates() const { return _ptree.subdivision_coordinates(); }

      /*! \brief The subdivision coordinates. */
      dimension_type subdivision_coordinate(size_type n) const { return _ptree.subdivision_coordinate(n); }

      /*! \brief The binary tree. */
      const BinaryTree& tree() const { return _ptree.tree(); }

      /*! \brief The mask. */
      const BooleanArray& mask() const { return _mask; }

      /*! \brief The number of cells in the set. */
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }

      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { 
        return const_iterator(scheme(),BinarySubtreeIterator(tree().begin(),_mask.begin())); 
      }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const {
        return const_iterator(scheme(),BinarySubtreeIterator(tree().end(),_mask.end())); 
      }

      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const {
        return tree().depth();
      }
      
      SizeArray subdivisions() const {
        SizeArray result(dimension(),1);
        for(size_type i=0; i!=depth(); ++i) {
          result[subdivision_coordinate(i)]*=2;
        }
        return result;
      }      
     private:
      PartitionTree<R> _ptree;
      BooleanArray _mask;
    };

  }
}

#endif /* _ARIADNE_PARTITION_TREE_SET_H */
