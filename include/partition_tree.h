/***************************************************************************
 *            partition_tree.h
 *
 *  1 July 2005
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file partition_grid.h
 *  \brief Cuboidal partition trees.
 */

#ifndef _PARTITION_TREE_H
#define _PARTITION_TREE_H

#include <algorithm>
#include <vector>
#include <iostream>

#include "array.h"
#include "sequence.h"
#include "utility.h"
#include "binary_word.h"

#include "rectangle.h"
#include "list_set.h"

namespace Ariadne {
  namespace Geometry {
      const bool leaf=1;
      const bool branch=0;
      const bool left=0;
      const bool right=1;

    typedef unsigned short dimension_type;
    typedef size_t size_type;
    typedef sequence<dimension_type> SubdivisionSequence;
    typedef std::vector<bool> BinaryArray;
    typedef std::vector<bool> BooleanArray;

    template<typename R> class Rectangle;
    template<typename R, template<typename> class BS> class ListSet;

    class BinaryTree;

    template<typename R> class PartitionScheme;
    template<typename R> class PartitionTreeCell;
    template<typename R> class PartitionTree;
    template<typename R> class PartitionTreeSet;

    template<typename R> class PartitionTreeIterator;
    template<typename R> class PartitionTreeSetIterator;

    class BinaryTreeIterator {
      friend class BinaryTree;
     private:
      BinaryTreeIterator(BinaryArray::const_iterator i) : _position(i), _word() { }
      BinaryTreeIterator& initialize() { while(*_position==branch) { _word.push_back(left); ++_position; } return *this; }
     public:
      BinaryTreeIterator(const BinaryArray& ba) : _position(ba.begin()), _word() { initialize(); }
     public:
      BinaryTreeIterator(const BinaryTreeIterator& iter) : _position(iter._position), _word(iter._word) { }
      bool operator==(const BinaryTreeIterator& iter) const { return this->_position==iter._position; }
      bool operator!=(const BinaryTreeIterator& iter) const { return this->_position!=iter._position; }
      const BinaryWord& operator*() const { return _word; }
      const BinaryWord* operator->() const { return &_word; }
      BinaryTreeIterator& operator++();
     private:
      BinaryArray::const_iterator _position;
      BinaryWord _word;
    };

    class BinaryTree {
     public:
      typedef BinaryTreeIterator const_iterator;
      BinaryTree() : _tree(1) { _tree[0]=leaf; }
      BinaryTree(const BinaryArray& t) : _tree(t) { }
      size_type capacity() const { return (_tree.size()+1)/2; }
      const_iterator begin() const { return const_iterator(_tree.begin()).initialize(); }
      const_iterator end() const { return const_iterator(_tree.end()); }
     private:
      BinaryArray _tree;
    };

    class BinaryTreeSetIterator {
     public:
      BinaryTreeSetIterator(const BinaryTree& t, const BooleanArray& m)
        : _tree_iter(t.begin()), _mask_iter(m.begin()) { initialize(); }
      BinaryTreeSetIterator(BinaryTree::const_iterator t,
                            BooleanArray::const_iterator m)
        : _tree_iter(t), _mask_iter(m) { }
      BinaryTreeSetIterator& initialize() { if(!*_mask_iter) { ++(*this); } return *this; }
      bool operator==(const BinaryTreeSetIterator iter) const {
        return _tree_iter==iter._tree_iter && _mask_iter==iter._mask_iter; }
      bool operator!=(const BinaryTreeSetIterator iter) const { return !(*this==iter); }
      const BinaryWord& operator*() const { return _tree_iter.operator*(); }
      const BinaryWord* operator->() const { return _tree_iter.operator->(); }
      BinaryTreeSetIterator& operator++() {
        do { ++_mask_iter; ++_tree_iter; } while((*_mask_iter==false) && !_tree_iter->empty()); return *this; }
     private:
      BinaryTreeIterator _tree_iter;
      BooleanArray::const_iterator _mask_iter;
    };

    /*!\brief A partition grid of rectangles in Euclidean space.
     */
    template<typename R> class PartitionScheme {
      typedef R real_type;
    public:
      /*! \brief Construct from a bounding box \a bb and a sequence of subdivision coordinates \a sc. */
      PartitionScheme(const Rectangle<R>& bb, const sequence<dimension_type>& sc)
        : _bounding_box(bb), _subdivision_coordinates(sc) { }

      bool operator==(const PartitionScheme<R>& pg) const {
        return _bounding_box==pg._bounding_box && _subdivision_coordinates==pg._subdivision_coordinates;
      }

      /*! \brief The underlying dimension of the grid. */
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


    /*! \brief A rectangle defined on a partition tree.
     */
    template<typename R>
    class PartitionTreeCell {
     public:
      typedef R real_type;
      typedef State<R> state_type;

      /*!\brief Construct from a grid and a binary word.
       */
      PartitionTreeCell(const PartitionScheme<R>& g,
                        const BinaryWord& w)
        : _scheme(g), _word(w)
      { }

      bool operator==(const PartitionTreeCell<R>& ptc) const {
        return _scheme==ptc._scheme && _word==ptc._word;
      }

      /*!\brief The underlying partition grid. */
      const PartitionScheme<R>& grid() const { return _scheme; }

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return _scheme.dimension(); }

      /*!\brief The dimension of the cell. */
      BinaryWord index() const { return _word; }

      /*!\brief Convert to an ordinary rectangle. */
      operator Rectangle<R>() const;
     private:
      PartitionScheme<R> _scheme;
      BinaryWord _word;
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

      /*! \brief Construct a one-cell tree. */
      PartitionTree(const PartitionScheme<R>& g)
        : _scheme(g), _tree() { }

      /*! \brief Construct a set based on a PartitionScheme, a subdivision tree and a mask. */
      PartitionTree(const PartitionScheme<R>& g, const BinaryArray& t)
        : _scheme(g), _tree(t) { }

      /*! \brief The underlying PartitionScheme. */
      const PartitionScheme<R>& scheme() const { return _scheme; }

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _scheme.dimension(); }

      /*! \brief The array describing the tree. */
      const BinaryTree& tree() const { return _tree; }

      /*! \brief The number of cells in the PartitionTree. */
      size_type capacity() const { return _tree.capacity(); }

      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const { return const_iterator(_scheme,_tree.begin()); }
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const { return const_iterator(_scheme,_tree.end()); }

     private:
      PartitionScheme<R> _scheme;
      BinaryTree _tree;
    };



   template<typename R>
    class PartitionTreeSetIterator {
      const PartitionScheme<R>& _scheme;
      BinaryTreeSetIterator _iter;
     public:
      PartitionTreeSetIterator(const PartitionScheme<R>& g, BinaryTreeSetIterator i) :
        _scheme(g), _iter(i) { }
      bool operator==(const PartitionTreeSetIterator<R>& other) { return this->_iter==other._iter; }
      bool operator!=(const PartitionTreeSetIterator<R>& other) { return this->_iter!=other._iter; }
      PartitionTreeCell<R> operator*() const { return PartitionTreeCell<R>(_scheme,*_iter); }
      PartitionTreeSetIterator& operator++() { ++_iter; return *this; }
      bool end() const { return _iter.end(); }
    };

  /*! \brief A denotable set on a partition grid, defined using a partition tree of cells.
     */
    template<typename R>
    class PartitionTreeSet {
     public:
      typedef PartitionTreeSetIterator<R> const_iterator;

      /*! \brief Construct an empty set based on a PartitionScheme. */
      PartitionTreeSet(const PartitionScheme<R>& g)
        : _ptree(g), _mask(1) { _mask[0]=false; }

      /*! \brief Construct a set based on a PartitionScheme, a subdivision tree and a mask. */
      PartitionTreeSet(const PartitionScheme<R>& g, const BinaryArray& t, const BooleanArray& m)
        : _ptree(g,t), _mask(m)
      {
        assert(_ptree.capacity() == _mask.size());
      }

      /*! \brief Construct a set based on a PartitionScheme, a subdivision tree and a mask. */
      PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m)
        : _ptree(t), _mask(m)
      {
        assert(_ptree.capacity() == _mask.size());
      }

      /*! \brief Convert to a ListSet of Rectangle. */
      operator ListSet<R,Rectangle> () const;

      /*! \brief The underlying PartitionScheme. */
      const PartitionScheme<R>& scheme() const { return _ptree.scheme(); }

      /*! \brief The space dimension of the set. */
      size_type dimension() const { return scheme().dimension(); }

      /*! \brief The binary tree. */
      const BinaryTree& tree() const { return _ptree.tree(); }

      /*! \brief The mask. */
      const BooleanArray& mask() const { return _mask; }

      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const {
        return const_iterator(scheme(),BinaryTreeSetIterator(tree().begin(),mask().begin()).initialize()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const {
        return const_iterator(scheme(),BinaryTreeSetIterator(tree().end(),mask().end())); }

     private:
      PartitionTree<R> _ptree;
      BooleanArray _mask;
    };




    template<typename R>
    PartitionTreeCell<R>::operator Rectangle<R>() const {
      Rectangle<R> res(_scheme.bounding_box());
      SubdivisionSequence::const_iterator coord_iter=_scheme.subdivision_coordinates().begin();
      BinaryWord::const_iterator word_iter=_word.begin();

      while(word_iter!=_word.end()) {
        size_type i=(*coord_iter);
        R centre = ( res.lower(i) + res.upper(i) ) / 2;
        if( (*word_iter)==left ) {
          res.set_upper(i,centre); }
        else {
          res.set_lower(i,centre);
        }
        ++word_iter;
        ++coord_iter;
      }

      return res;
    }


    template<typename R>
    PartitionTreeSet<R>::operator ListSet<R,Rectangle>() const {
      ListSet<R,Rectangle> res(this->dimension());
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(Rectangle<R>(*iter));
      }
      return res;
    }


    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionScheme<R>& g)
    {
      os << "{ class: PartitionScheme<" << name<R>() << ">,\n";
      os << "  bounding_box: " << bounding_box() << ",\n";
      os << "  subdivision_coordinates: " << subdivision_coordinates() << "\n";
      os << "}\n";
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeCell<R>& c)
    {
      os << "{ class: PartitionTreeCell<" << name<R>() << ">,\n";
      os << "  word: " << c.index() << ",\n";
      os << "  rectangle: " << Rectangle<R>(c) << "\n";
      os << "}\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTree<R>& t)
    {
      os << "{ class: PartitionTree<" << name<R>() << ">,\n";
      os << "  words: ";
      write_array(os, t.tree().begin(), t.tree().end());
      os << "\n";
      os << "  rectangles: [ " << Rectangle<R>(*t.begin());
      for(typename PartitionTree<R>::const_iterator ptree_iter=++t.begin() ; ptree_iter!=t.end(); ++ptree_iter) {
        os << ", " << Rectangle<R>(*ptree_iter);
      }
      os << " ]\n";
      os << "}\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeSet<R>& set)
    {
      os << "{ class: PartitionTreeSet<" << name<R>() << ">,\n";
      os << "  bounding_box: " << set.scheme().bounding_box() << ",\n";
      os << "  subdivision_coordinates: " << set.scheme().subdivision_coordinates() << ",\n";
      os << "  words: ";
      write_array(os, BinaryTreeSetIterator(set.tree(),set.mask()), BinaryTreeSetIterator(set.tree().end(),set.mask().end()));
      os << "\n";
      os << "  rectangles: [ " << Rectangle<R>(*set.begin());
      for(typename PartitionTreeSet<R>::const_iterator pts_iter=++set.begin() ; pts_iter!=set.end(); ++pts_iter) {
        os << ", " << Rectangle<R>(*pts_iter);
      }
      os << " ]\n";
      os << "}\n";
      return os;
    }

  }
}

#include "partition_tree.cc"

#endif /* _PARTITION_TREE_H */
