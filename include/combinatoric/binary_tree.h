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

#ifndef ARIADNE_BINARY_TREE_H
#define ARIADNE_BINARY_TREE_H

#include <vector>
#include <iosfwd>

#include <boost/iterator/iterator_facade.hpp>

#include "../base/types.h"
#include "../base/array.h"

#include "../combinatoric/binary_word.h"

namespace Ariadne {
  namespace Combinatoric {
    class BinaryTree;
    class BinaryTreeIterator;
    
    /*!\brief A sorted list of BinaryWord elements such that every binary sequence has a unique prefix in the list.
     *   Optimised for memory usage.
     *   Constant forward iterators only, access may be inefficient.
     */
    class BinaryTree {
     public:
      /*!\brief The constant indicating that a node is a leaf of the tree. */
      static const bool leaf=1;
      /*!\brief The constant indicating that a node is a branch of the tree. */
      static const bool branch=0;
      /*!\brief The constant indicating that a word follows the left-hand path. */
      static const bool left=0;
      /*!\brief The constant indicating that a word follows the right-hand path. */
      static const bool right=1;
     public:    
      /*! \brief The type of the BinaryWord object stored in the tree. */
      typedef BinaryWord value_type;
      typedef BinaryTreeIterator iterator;
      typedef BinaryTreeIterator const_iterator;
      friend class BinaryTreeIterator;
     public:
      /*!\brief Construct a tree containing just the empty word. */
      explicit BinaryTree();
      /*!\brief Construct a tree of depth n containing all words of size n. */
      explicit BinaryTree(size_type n);

      /*! \brief Construct a tree from an array of bits denoting branches and leaves. */
      explicit BinaryTree(const array<bool>& t);
      explicit BinaryTree(const std::vector<bool>& t);
      explicit BinaryTree(const size_type& n, const bool* ptr);
      
      /*! \brief The total number of leaf nodes in the tree. */
      size_type size() const;
      
      /*! \brief A constant forward iterator to the first leaf of the tree. */
      const_iterator begin() const;
      /*! \brief A constant forward iterator to the end of the tree. */
      const_iterator end() const;
      
      /*! \brief The depth of the furthest leaf of the tree. */
      size_type depth() const;
      
      /*! \brief The array of bits denoting branches and leaves. */
      const std::vector<bool>& array() const;
      
      /*!\brief Advance a binary word used to describe a position in a tree by
       * removing all trailing right's and replacing the last left by a right.    
       */
      static void advance(BinaryWord& word);
      
      friend std::ostream& operator<<(std::ostream&, const BinaryTree&);
     private:
      void check() const;
     private:
      std::vector<bool> _array;
    };
    
    class BinaryTreeIterator 
      : public boost::iterator_facade<BinaryTreeIterator,
                                      BinaryWord,
                                      boost::forward_traversal_tag,
                                      const BinaryWord&>
    {
      friend class BinaryTree;
     private:
      BinaryTreeIterator(BooleanArray::const_iterator i, BooleanArray::const_iterator e);
     public:
      void skip_subtree();
     private:
      bool equal(const BinaryTreeIterator& other) const { 
        return this->_position==other._position && this->_word==other._word; }
      const BinaryWord& dereference() const { 
        return _word; }
      void increment();
      friend class boost::iterator_core_access;
     private:
      // Invariant: _position should always point to a leaf after operator++
      std::vector<bool>::const_iterator _position;
      BinaryWord _word;
    };
    

    inline
    size_type 
    BinaryTree::size() const
    { 
      return (_array.size()+1)/2; 
    }
      
    inline
    const std::vector<bool>& 
    BinaryTree::array() const 
    { 
      return this->_array;
    }
    
    inline
    BinaryTree::const_iterator 
    BinaryTree::begin() const 
    { 
      return const_iterator(_array.begin(),_array.end());
    }
    
    inline
    BinaryTree::const_iterator 
    BinaryTree::end() const 
    { 
      return const_iterator(_array.end(),_array.end()); 
    }
    
    std::ostream& operator<<(std::ostream& os, const BinaryTree& t);
    
  } // namespace Base
}// namespace Ariadne


#endif /* ARIADNE_BINARY_TREE_H */
