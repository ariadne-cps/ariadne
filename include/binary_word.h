/***************************************************************************
 *            binary_word.h
 *
 *  18 January 2005
 *  Copyright  2004,2005  Alberto Casagrande, Pieter Collins
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
 
/*! \file binary_word.h
 *  \brief Binary words and sets of binary words.
 *
 * Binary words are generally useful objects. 
 * This file provides a class of binary words of 
 * arbitrary size, complementing the STL bitset and
 * vector<bool> classes. The main purpose of giving 
 * this class is to have a known implementation of 
 * binary words which may help in giving efficient
 * implementations of the list classes described below.
 * The most efficient implementation may, in fact, be an
 * STL vector<bool>. Hence another purpose of this class
 * is to provide a known interface for binary words.
 *
 * We give a number of classes representing sets of binary words. 
 * These classes have been designed to be optimised for storage, and for 
 * sequential iteration through the elements. Hence the preferred method of
 * access is through iterators, which need only conform to the requirements
 * of an InputIterator. These  classes may be immutable. 
 * The fundamental operations are to join two sets,
 * and to convert between different representations.
 * 
 * A special class is given representing prefix-free sets of words
 * using a binary tree representation. 
 *
 * The main intended use of these classes in Ariadne is to represent 
 * sets as unions of rectangles on a grid. See the file grid_rectangle.h
 * for the implementation of this representation. 
 * The basic idea is that a subset of Euclidean space contained in some 
 * cuboid R can be represented by a union of sets obtained by repeatedly
 * subdividing R in two along the coordinate axes. 
 * The nth subdivision occurs along a coordinate specified by the Grid 
 * class, and whether the lower or upper set is taken depends the value 
 * of the nth element of a BinaryWord. 
 * In this way, sets of points can be represented efficiently by just listing
 * in which half of the remaining rectangle the subdivision occurs in.
 *
 * NOTE TO ALBERTO : This describes the basic ideas. Hope it's comprehensible. 
 * I haven't finished the implementation of all features. There are loose ends
 * in how we choose subdivision directions. I hope all this effort will result
 * in a better implementation than "pointer-based" quad trees; it should be
 * much more memory-efficient, but may be unacceptably slow. On the other hand,
 * for big grids or fine subdivisions, memory may well be a limiting factor.
 */

#ifndef _BINARY_WORD_H
#define _BINARY_WORD_H

#include <bitset>
#include <vector>
#include <exception>
#include <iostream>

#include "utility.h"

namespace Ariadne 
{	
  class _BinaryWordList_const_iterator;
  class _BinaryWordixedSizeList_const_iterator;


  /*! \brief A statically-allocated binary word of fixed maximum length.
   */
  class BinaryWord {
    friend class BinaryWordList;
    friend class BinaryWordPrefixTree;
    friend class BinaryWordPrefixFreeTree;
    friend class BinaryWordFixedSizeList;
   public:
    /*! \brief An unsigned integral type. */
    typedef std::vector<bool>::size_type size_type;
    typedef std::vector<bool>::const_iterator const_iterator;
   private:
    /* The type of a machine byte (not used) */
    typedef unsigned char byte_type; 
    /* The number of bits per byte */
    static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
   private:
    std::vector<bool> _rep;
   public:
    /*! \brief Default constructor creates a zero-length word.  */
    BinaryWord() : _rep() { }
    
    /*! \brief Conversion from a std:vector<bool>. */
    BinaryWord(const std::vector<bool>& v) : _rep(v) { }
    
    /*! \brief Create a word of length \a n from the boolean array starting at \a p. */
    BinaryWord(size_type n, const bool* p) : _rep(n,p) { }
    
    /*! \brief Create a word of length \a n from the boolean array starting at \a iter. */
    BinaryWord(size_type n, std::vector<bool>::const_iterator first) : _rep(first,first+n) { }
    
    /*! \brief Create a word of length \a n from the boolean array between \a first and \a last. */
    BinaryWord(std::vector<bool>::const_iterator first, std::vector<bool>::const_iterator last) : _rep(first,last) { }
    
    /*! \brief Create a word of length \a n from the chunk of memory starting at \a p. */
    BinaryWord(size_type n, const void* p) : _rep(n,reinterpret_cast<const bool*>(p)) { }
    
    /*! \brief Is the word empty. */
    bool empty() const { return _rep.empty(); }

    /*! \brief The size of the word, in bits. */
    size_type size() const { return _rep.size(); }

    /*! \brief The largest possible size of the BinaryWord. */
    size_type max_size() const { return _rep.max_size(); }
    
    /* \brief The number of machine bytes the word takes up. */
    size_type bytes() const { return (size()+_bits_per_byte-1)/_bits_per_byte; }
    
    /*! \brief Returns the nth bit. */
    bool operator[] (size_type n) const { return _rep[n]; }

    /*! \brief Checked access to the nth bit. */
    bool at(size_type n) const { if(n<size()) { return operator[](n); } else { throw std::out_of_range("BinaryWord"); } }

    /*! \brief Returns the last bit. */
    bool back () const { return operator[](size()-1); }

     /*! \brief Sets the last bit to \a x.  */
    void set_back (bool x) { set(size()-1,x); }
 
    /*! \brief Insert a bit at the end of the word. */
    void push_back(bool x) { _rep.push_back(x); }
    
    /*! \brief Removes the bit at the end of the word. */
    void pop_back() { _rep.pop_back(); }
    
    /* \brief An iterator to the beginning of the bits of the BinaryWord. */
    const_iterator begin() const { return _rep.begin(); }
    
    /* \brief An iterator to the beginning of the bits of the BinaryWord. */
    const_iterator end() const { return _rep.end(); }
    
    /*! \brief true if the word is a subword of the other word. */
    bool is_subword(const BinaryWord& b) const {
      if(size()>b.size()) { return false; }
      for(size_type i=0; i!=size(); ++i) { if((*this)[i] != b[i]) { return false; } }
      return true;
    }

    /*! \brief Stream insertion operator. */
    friend std::ostream& operator<<(std::ostream& os, const BinaryWord& bw);
    /*! \brief Stream extraction operator. */
    friend std::istream& operator>>(std::istream& is, BinaryWord& bw);
   private:
    /* Set the nth element to x. */
    void set(const size_type& n, const bool& x) { _rep[n]=x; }
  };


  /*!\brief A list of BinaryWord elements of variable size, optimised for memory usage. Access by constant forward iterators only.
   */
  class BinaryWordList {
    /*!\brief Union of two BinaryWordList s. */
    friend BinaryWordList join(const BinaryWordList& l1, const BinaryWordList& l2);
    
    typedef size_t size_type;
    typedef BinaryWord::byte_type byte_type;
    typedef _BinaryWordList_const_iterator const_iterator;
   private:
    static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
    std::vector<bool> _elements;
    std::vector<size_type> _sizes;
   public:
    /*!\brief Construct an empty list.
     */
    BinaryWordList() : _sizes(), _elements() { }
    
    /*!\brief The size of the list, in binary words. */
    size_type size() const { return _sizes.size(); }
    
    /*!\brief Insert an element at the back of the list. */
    void push_back(const BinaryWord& b) { 
      _sizes.push_back(b.size());
      for(int i=0; i!=b.size(); ++i) { _elements.push_back(b[i]); }
    }
    
    /*!\brief Returns a constant forward iterator pointing to the beginning of the BinaryWordList. */
    const_iterator begin() const;

    /*!\brief Returns a constant forward iterator pointing to the end of the BinaryWordList. */
    const_iterator end() const;
  };
  
  /*!\brief A list of BinaryWord elements, each of the same size, optimised for memory usage. */
  class BinaryWordFixedSizeList {
    friend class _BinaryWordFixedSizeList_const_iterator;
    
    /*!\brief Union of two BinaryWordFixedSizeList s.
     */
    friend BinaryWordFixedSizeList join(const BinaryWordFixedSizeList& l1, const BinaryWordFixedSizeList& l2);
    
    typedef size_t size_type;
    typedef BinaryWord::byte_type byte_type;
    typedef  _BinaryWordFixedSizeList_const_iterator const_iterator;
   private:
    static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
    size_type _word_size;
    std::vector<bool> _elements;
   public:
    BinaryWordFixedSizeList(size_type ws) : _word_size(ws), _elements() { }
    
    /*!\brief The size of the list, in binary words.
     */
    size_type size() const { return _elements.size()/word_size(); }

    /*!\brief The size of each binary word in the list, in bits. */
    size_type word_size() const { return _word_size; }
    
    /*!\brief Insert an element at the back of the list. */
    void push_back(const BinaryWord& b) { 
      assert(b.size()==word_size()); 
      for(int i=0; i!=word_size(); ++i) { _elements.push_back(b[i]); }
    }

    /*!\brief Remove the element at the back of the list. */
    void pop_back() { for(int i=0; i!=word_size(); ++i) { _elements.pop_back(); } }

    /*!\brief Returns the element at the back of the list. */
    BinaryWord back() const {
      std::vector<bool>::const_iterator iter=_elements.end()-word_size();
      return BinaryWord(word_size(), iter);
    }
    
    /*!\brief Returns the nth element of the list.
     */
    BinaryWord operator[] ( size_type n) const {
      std::vector<bool>::const_iterator iter=_elements.begin() + (n*word_size());
      return BinaryWord(word_size(), iter);
    }
    
    /*!\brief Returns a constant random-access iterator pointing to the beginning of the BinaryWordFixedSizeList.  */
    const_iterator begin() const;

    /*!\brief Returns a constant random-access iterator pointing to the end of the BinaryWordFixedSizeList. */
    const_iterator end() const;
    
   private:
    /* \brief The number of bytes needed to store the list. */
    size_type storage_in_bytes() const { return _elements.size()/_bits_per_byte; }
  };
  
  /*!\brief A sorted list of BinaryWord elements such that every binary sequence has a unique prefix in the list.
   *   Optimised for memory usage.
   *   Constant forward iterators only, access may be inefficient.
   */
  class BinaryWordPrefixTree {
    friend class _BinaryWordPrefixTree_const_iterator;

    typedef size_t size_type;
    typedef BinaryWord::byte_type byte_type;
   private:
    std::vector<bool> _rep;
    
    static const bool _leaf = 1;
    static const bool _branch = 0;
   public:
    /*! \brief The type of the BinaryWord object stored in the tree. */
    typedef BinaryWord value_type;

    typedef _BinaryWordPrefixTree_const_iterator const_iterator;
    
   public:
    /*!\brief Construct a list containing just the empty word.
     */
    BinaryWordPrefixTree() : _rep() { _rep.push_back(_leaf); }
    
    BinaryWordPrefixTree(const std::vector<bool>& __rep) : _rep(__rep) { }
    BinaryWordPrefixTree(const size_type& __n, const bool* __ptr) : _rep(__n,__ptr) { }
    
    /*!\brief The number of words in the list.
     */
    size_type size() const { return (_rep.size()+1)/2; }
    
    /*!\brief Returns a constant forward iterator pointing to the beginning of the BinaryWordPrefixTree.
     */
    const_iterator begin() const;
    /*!\brief Returns a constant forward iterator pointing to the end of the BinaryWordPrefixTree.
     */
    const_iterator end() const;
  };
    
  class _BinaryWordPrefixTree_const_iterator {
    typedef _BinaryWordPrefixTree_const_iterator Self;
    typedef size_t size_type;
    
    std::vector<bool>::const_iterator _pos;
    BinaryWord _word;
      
    static const bool _leaf = BinaryWordPrefixTree::_leaf;
    static const bool _branch = BinaryWordPrefixTree::_branch;
   public:
    _BinaryWordPrefixTree_const_iterator(const std::vector<bool>::const_iterator& __pos)
      : _pos(__pos), _word() { }
      
    _BinaryWordPrefixTree_const_iterator(const Self& __other) 
      : _pos(__other._pos), _word(__other._word) { }
      
    bool operator==(const Self& __other) const {
      return (_pos==__other._pos); 
    }
      
    BinaryWord operator*() const { _word; }
      
    Self& operator++() { 
      if( !_word.empty() ) {
	while( _word.back() == 1 ) {
	  _word.pop_back();
	  if( _word.empty() ) {
	    return *this;
	  }
	}
	_word.set_back(1);
      }
      while( (*_pos) == _branch ) {
	++_pos;
	_word.push_back(0);
      }
      
      ++_pos;
      return *this;
    }
      
    Self operator++(int) {
      Self res=(*this);
      ++(*this);
      return res;
    }
  };

  inline BinaryWordPrefixTree::const_iterator BinaryWordPrefixTree::begin() const { return ++(const_iterator(_rep.begin())); }
  inline BinaryWordPrefixTree::const_iterator BinaryWordPrefixTree::end() const { return const_iterator(_rep.end()); }
      
  
  /*!\brief A list of prefix-free BinaryWord elements using a tree representation.
   *   Optimised for memory usage.
   *   Constant forward iterators only, access may be inefficient.
   */
  class BinaryWordPrefixFreeTree {
    friend class _BinaryWordPrefixFreeTree_const_iterator;
    /*!\brief Union of two BinaryWordPrefixFreeTree s.
     */
    friend BinaryWordPrefixFreeTree join(const BinaryWordPrefixFreeTree& __tree1, const BinaryWordPrefixFreeTree& __tree2);
    
    typedef size_t size_type;
    typedef BinaryWord::byte_type byte_type;
   private:
    static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
    BinaryWordPrefixTree _tree;
    std::vector<bool> _mask;
    
    static const bool _leaf = 0;
    static const bool _branch = 1;
    static const bool _empty = 0;
    static const bool _full = 1;
   public:
    /* \brief A forward iterator through the tree.
     */
    typedef _BinaryWordPrefixFreeTree_const_iterator const_iterator;
    
    /*!\brief Construct an empty list of prefix-free binary words 
     */
    BinaryWordPrefixFreeTree() : _tree(), _mask() { }
    /*    
    BinaryWordPrefixFreeTree(const size_type& __n, const std::vector<bool>& __rep) : _size(__n), _rep(__rep) { }
    BinaryWordPrefixFreeTree(const size_type& __n, const bool* __ptr) : _size(__n), _rep(__n,__ptr) { }
    BinaryWordPrefixFreeTree(const BinaryWordList& __list);
    */

    /*!\brief The number of words in the list.
     */
    size_type size() const { return std::count(_mask.begin(), _mask.end(), true); }
    
    /*!\brief Returns a constant forward iterator pointing to the beginning of the BinaryWordPrefixFreeTree.
     */
    const_iterator begin() const;
    /*!\brief Returns a constant forward iterator pointing to the end of the BinaryWordPrefixFreeTree.
     */
    const_iterator end() const;
   private:
    /* \brief The number of bytes needed to store the list.
     */
    // size_type storage_in_bytes() const { return _rep.size() / _bits_per_byte; }
  };
  
  class _BinaryWordPrefixFreeTree_const_iterator {
    typedef _BinaryWordPrefixFreeTree_const_iterator Self;
    typedef size_t size_type;
    
    std::vector<bool>::const_iterator _pos;
    std::vector<bool> _word;
    
    static const bool _leaf = BinaryWordPrefixFreeTree::_leaf;
    static const bool _branch = 1;
    static const bool _empty = 0;
    static const bool _full = 1;
   public:
    _BinaryWordPrefixFreeTree_const_iterator(const std::vector<bool>::const_iterator& __pos)
      : _pos(__pos), _word() {
    }
    
    _BinaryWordPrefixFreeTree_const_iterator(const Self& __other) 
      : _pos(__other._pos), _word(__other._word) {
    }
    
    bool operator==(const Self& __other) const {
      return (_pos==__other._pos); 
    }
    
    BinaryWord operator*() const { return BinaryWord(_word); }
    
    Self& operator++() { 
      while(true) {
	if( !_word.empty() ) {
	  while( _word.back() == 1 ) {
	    _word.pop_back();
	    if( _word.empty() ) {
	      return *this;
	    }
	  }
	  _word.back()=1;
	}
	while( (*_pos) == _branch ) {
	  ++_pos;
	  _word.push_back(0);
	}
	if( (*_pos) == _full ) {
	  ++_pos;
	  return *this;
	} 
	else {
	  ++_pos;
	}
      }
    }
    
    Self operator++(int) {
      Self res=(*this);
      ++(*this);
      return res;
    }
    
    void skip_subtree() {
      size_type __len=_word.size();
      do {
	if(*_pos == _branch) {
	  _word.push_back(0);
	  ++_pos;
	}
	else {
	  ++_pos;
	  ++_pos;
	  if(_word.size()==__len) {
	    return;
	  }
	  while(_word.back()==1) {
	    _word.pop_back();
	    if(_word.size()==__len) {
	      return;
	    }
	  }
	  _word.pop_back();
	  _word.push_back(1);
	}
      }
      while(true);
    }
  };

  /*
  
  inline BinaryWordPrefixFreeTree::const_iterator BinaryWordPrefixFreeTree::begin() const { return ++(const_iterator(_rep.begin())); }
  inline BinaryWordPrefixFreeTree::const_iterator BinaryWordPrefixFreeTree::end() const { return const_iterator(_rep.end()); }

  */
  
  BinaryWordPrefixFreeTree 
  inline
  join(const BinaryWordPrefixFreeTree& __tree1, const BinaryWordPrefixFreeTree& __tree2) {
    BinaryWordPrefixFreeTree result;
    BinaryWordPrefixTree& _tree = result._tree;
    std::vector<bool>& mask = result._mask;

    BinaryWord word;
    /*
      BinaryWordPrefixFreeTree::const_iterator iter1=__tree1.begin();
      BinaryWordPrefixFreeTree::const_iterator iter2=__tree2.begin();
    */	
    BinaryWordPrefixFreeTree::const_iterator iter1=__tree1.begin();
    BinaryWordPrefixFreeTree::const_iterator iter2=__tree2.begin();
    
    static const bool _leaf = BinaryWordPrefixFreeTree::_leaf;
    static const bool _branch = BinaryWordPrefixFreeTree::_branch;
    static const bool _empty = BinaryWordPrefixFreeTree::_empty;
    static const bool _full =  BinaryWordPrefixFreeTree::_full;
    
   
    
    /* FIXME!

    bool bool1=*iter1;
    bool bool2=*iter2;
    
    if( ((*iter1)==_branch) or ((*iter2)==_branch) ) {
      tree.push_back(_branch); 
      ++iter1;
      ++iter2;
    }
    else if ( (*iter1==_leaf) and (*iter2==_leaf) ) {
      tree.push_back(_leaf); 
      ++iter1;
      ++iter2;
      if( (*iter1==_full) or (*iter2==_full) ) {
	tree.push_back(_full);
      }
      else {
	tree.push_back(_empty);
      }
      ++iter1;
      ++iter2;
    }
    else if( (*iter1==_branch) and (*iter2==_branch) ) {
      tree.push_back(_branch);
      ++iter1;
      ++iter2;
    }
    else {
      if( *iter1==_branch ) {
	++iter2;
	// BinaryWordPrefixFreeTree::const_iterator iter0=iter1;
	std::vector<bool>::const_iterator iter0=iter1;
	// FIXME -- this is for tree iterators!
	// iter1.skip_subtree
	if(*iter2==_full) {
	  tree.push_back(_leaf);
	  tree.push_back(_full);
	}
	else {
	  while(iter0!=iter1) {
	    tree.push_back(*iter0);
	    ++iter0;
	  }
	}
      }
      else { // iter2==_branch
	++iter1;
	// BinaryWordPrefixFreeTree::const_iterator iter0=iter2;
	std::vector<bool>::const_iterator iter0=iter2;
	// FIXME
	// iter2.skip_tree;
	if(*iter1==_full) {
	  tree.push_back(_leaf);
	  tree.push_back(_full);
	}
	else {
	  while(iter0!=iter2) {
	    tree.push_back(*iter0);
	    ++iter0;
	  }
	}
      }
    }
    if( (tree[tree.size()-3] == _branch) && (tree[tree.size()-2]==tree[tree.size()-1]) ) {
      bool state=tree.back();
      tree.pop_back();
      tree.pop_back();
      tree.pop_back();
      tree.push_back(state);
    }
    */
  }
  
  
  inline
  std::istream& operator>>(std::istream& is, BinaryWord& b) 
  {    
    std::vector<bool> v;
    is >> v;
    b=BinaryWord(v);

    return is;
  }
  
  inline
  std::ostream& operator<<(std::ostream& os, const BinaryWord& b) 
  {    
    os << b.size() << " :"; for(int i=0; i!=b.size(); ++i) { if(i%8==0) { os << " "; } os << b[i];  }
    return os;
  }
  



}

#endif /* _BINARY_WORD_H */

