/***************************************************************************
 *            binary_tree.h
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

/*! \file binary_tree.h
 *  \brief Cuboidal binary trees.
 */

#ifndef _ARIADNE_BINARY_TREE_H
#define _ARIADNE_BINARY_TREE_H

#include <algorithm>
#include <vector>
#include <iostream>

#include "basic_type.h"
#include "binary_word.h"

#include "utility.h"

namespace Ariadne {
  const bool leaf=1;
  const bool branch=0;
  const bool left=0;
  const bool right=1;
  
  class BinaryTree;
  class BinaryTreeIterator;

  class BinaryTreeIterator {
    friend class BinaryTree;
   private:
    BinaryTreeIterator(BooleanArray::const_iterator i) : _position(i), _word() { }
    BinaryTreeIterator& initialize() { while(*_position==branch) { _word.push_back(left); ++_position; } return *this; }
   public:
    BinaryTreeIterator(const BooleanArray& ba) : _position(ba.begin()), _word() { initialize(); }
   public:
    BinaryTreeIterator(const BinaryTreeIterator& iter) : _position(iter._position), _word(iter._word) { }
    bool operator==(const BinaryTreeIterator& iter) const { return this->_position==iter._position; }
    bool operator!=(const BinaryTreeIterator& iter) const { return this->_position!=iter._position; }
    const BinaryWord& operator*() const { return _word; }
    const BinaryWord* operator->() const { return &_word; }
    BinaryTreeIterator& operator++();
    void skip_subtree();
   private:
    // Invariant: _position should always point to a leaf after operator++
    BooleanArray::const_iterator _position;
    BinaryWord _word;
  };

  /*!\brief A sorted list of BinaryWord elements such that every binary sequence has a unique prefix in the list.
   *   Optimised for memory usage.
   *   Constant forward iterators only, access may be inefficient.
   */
  class BinaryTree {
    friend class BinaryTreeIterator;
  
    typedef BinaryWord::byte_type byte_type;
   public:
    /*! \brief The type of the BinaryWord object stored in the tree. */
    typedef BinaryWord value_type;

    typedef size_t size_type;
    typedef BinaryTreeIterator const_iterator;
   public:
    /*!\brief Construct a list containing just the empty word. */
    explicit BinaryTree() : _array(1) { _array[0]=leaf; }
    /*! \brief Construct a tree from an array of bits denoting branches and leaves. */
    explicit BinaryTree(const BooleanArray& t) : _array(t) { check(); }
    explicit BinaryTree(const size_type& n, const bool* ptr) : _array(n,ptr) { check(); }
  
    /*! \brief The total number of leaf nodes in the tree. */
    size_type size() const { return (_array.size()+1)/2; }
    /*! \brief A constant forward iterator to the first leaf of the tree. */
    const_iterator begin() const { return const_iterator(_array.begin()).initialize(); }
    /*! \brief A constant forward iterator to the end of the tree. */
    const_iterator end() const { return const_iterator(_array.end()); }
    
    friend std::ostream& operator<<(std::ostream&, const BinaryTree&);
   private:
    void check() const { assert(2u*std::count(_array.begin(),_array.end(),branch)+1u==_array.size()); }
   private:
    BooleanArray _array;
  };
    
  std::ostream& operator<<(std::ostream& os, const BinaryTree& t);
 
}

  
#endif /* _ARIADNE_BINARY_TREE_H */
