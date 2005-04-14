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

namespace Ariadne 
{	
    class _BinaryWord_const_iterator;

    /*! \brief A statically-allocated binary word of fixed maximum length.
     */
    class BinaryWord {
	friend class BinaryWordList;
	friend class BinaryWordPrefixFreeTree;
	friend class BinaryWordFixedSizeList;
      public:
	/*! \brief An unsigned integral type.
	 */
	typedef size_t size_type;
	
	/*! \brief  The type of byte used to store the word.
	 */
	typedef unsigned char byte_type;

	/*! \brief  The type of byte used to store the word.
	 */
	typedef _BinaryWord_const_iterator const_iterator;
      private:
	/* The number of bits per byte */
	static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
	static const size_type _storage_size_in_bytes=32;
	byte_type _rep[_storage_size_in_bytes];
      public:
	/*! \brief Default constructor creates a zero-length word. 
	 */
	BinaryWord() { _rep[0]=0; };
	
	/*! \brief Conversion from a std:vector<bool>.
	 */
	BinaryWord(const std::vector<bool>& v);

	/*! \brief Create a word of length \a n from the chunk of memory starting at \a p.
	 */
	BinaryWord(size_type n, const void* p) { 
	    if(n>max_size()) { throw std::invalid_argument("Binary word size"); } 
	    _rep[0]=n; const byte_type* q=(byte_type*) p; for(size_type i=0; i!=bytes(); ++i) { _rep[i+1]=q[i]; } 
	}
	
	//    BinaryWord(size_type n, const std::bitset<256>& bs) : _size(n), _bits(bs) { }
	
	/*! \brief The size of the word, in bits.
	 */
	size_type size() const { return size_type(_rep[0]); }
	/*! \brief The largest possible size of the BinaryWord.
	 */
	size_type max_size() const { return (_storage_size_in_bytes-1)*_bits_per_byte; }
	
	/* \brief The number of machine bytes the word takes up. */
	size_type bytes() const { return (size()+_bits_per_byte-1)/_bits_per_byte; }
	
	/*! \brief Returns the nth bit.
	 */
	bool operator[] (size_type n) const { byte_type b=_rep[n/_bits_per_byte+1]; return bool(b & (1<<(n%_bits_per_byte))); }
	/*! \brief Checked access to the nth bit.
	 */
	bool at(size_type n) const { if(n<size()) { return operator[](n); } else { throw std::out_of_range("BinaryWord"); } }
	
	/* \brief A pointer to the beginning of the bits of the BinaryWord. FIXME: make private! 
	 */
	const_iterator begin() const;
	
	/* \brief A pointer to the beginning of the bits of the BinaryWord. FIXME: make private! 
	 */
	const void* pointer() const { return ((_rep)+1); }    
	
	/* \brief A pointer to the beginning of the bits of the BinaryWord. FIXME: make private! 
	 */
	const void* pointer() const { return ((_rep)+1); }    
	
	/*! \brief true if the word is a subword of the other word. 
	 */
	bool is_subword(const BinaryWord& b) const {
	    if(size()>b.size()) { return false; }
	    for(size_type i=0; i!=size(); ++i) { if((*this)[i] != b[i]) { return false; } }
	    return true;
	}
    };
    
    /*!\brief A list of BinaryWord elements of variable size, optimised for memory usage. Access by constant forward iterators only.
     */
    class BinaryWordList {
	friend class _BinaryWordList_const_iterator;
	/*!\brief Union of two BinaryWordList s.
	 */
	friend BinaryWordList join(const BinaryWordList& __list1, const BinaryWordList& __list2);

	typedef size_t size_type;
	typedef BinaryWord::byte_type byte_type;
	typedef _BinaryWordList_const_iterator const_iterator;
      private:
	static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
	size_type _size;
	std::vector<byte_type> _rep;
      public:
	/*!\brief Construct an empty list.
	 */
	BinaryWordList() : _size(0), _rep() { }
	
	/*!\brief The size of the list, in binary words.
	 */
	size_type size() const { return _size; }
	
	/*!\brief Insert an element at the back of the list.
	 */
	void push_back(const BinaryWord& b) { 
	    ++_size; _rep.push_back( (byte_type) (b.size()) ); 
	    for(int i=0; i!=b.bytes(); ++i) { _rep.push_back(((byte_type*) b.pointer()) [i] ); } 
	}

	/*!\brief Returns a constant forward iterator pointing to the beginning of the BinaryWordList.
	 */
	const_iterator begin() const;
	/*!\brief Returns a constant forward iterator pointing to the end of the BinaryWordList.
	 */
	const_iterator end() const;
      private:
	/* \brief The number of bytes needed to store the list.
	 */
	size_type storage_in_bytes() const { return _rep.size(); }
    };
    
    /*!\brief A list of BinaryWord elements, each of the same size, optimised for memory usage.
     */
    class BinaryWordFixedSizeList {
	friend class _BinaryWordFixedSizeList_const_iterator;

	/*!\brief Union of two BinaryWordFixedSizeList s.
	 */
	friend BinaryWordFixedSizeList join(const BinaryWordFixedSizeList& __list1, const BinaryWordFixedSizeList& __list2);

	typedef  size_t size_type;
	typedef  BinaryWord::byte_type byte_type;
      private:
	static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
	size_type _word_size;
	std::vector<byte_type> _rep;
	size_type _word_bytes() const { return (_word_size+_bits_per_byte-1)/_bits_per_byte; }
      public:
	BinaryWordFixedSizeList(size_type ws) : _word_size(ws), _rep() { }
	
	/*!\brief The size of the list, in binary words.
	 */
	size_type size() const { return _rep.size()/_word_bytes(); }
	/*!\brief The size of each binary word in the list, in bits.
	 */
	size_type word_size() const { return _word_size; }
	
	/*!\brief Insert an element at the back of the list.
	 */
	void push_back(const BinaryWord& b) { 
	    assert(b.size()==word_size()); 
	    for(int i=0; i!=_word_bytes(); ++i) { 
		_rep.push_back(((byte_type*) b.pointer()) [i] ); } 
	}
	/*!\brief Remove the element at the back of the list.
	 */
	void pop_back() { 
	    for(int i=0; i!=_word_bytes(); ++i) { _rep.pop_back(); } }
	/*!\brief Returns the element at the back of the list.
	 */
	BinaryWord back() const {
	    const byte_type* __ptr=&_rep[(size()-1)*_word_bytes()];
	    return BinaryWord(_word_size, __ptr);
	}

	/*!\brief Returns the nth element of the list.
	 */
	BinaryWord operator[] ( size_type __n) const {
	    const byte_type* __ptr=&_rep[__n*_word_bytes()];
	    return BinaryWord(_word_size, __ptr);
	}

	/*!\brief Returns a constant random-access iterator pointing to the beginning of the BinaryWordFixedSizeList.
	 */
	const_iterator begin() const;
	/*!\brief Returns a constant random-access iterator pointing to the end of the BinaryWordFixedSizeList.
	 */
	const_iterator end() const;

      private:
	/* \brief The number of bytes needed to store the list.
	 */
	size_type storage_in_bytes() const { return _rep.size(); }
    };
    
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
	size_type _size;
	std::vector<bool> _rep;
 
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
	BinaryWordPrefixFreeTree() : _size(0), _rep() { }

	BinaryWordPrefixFreeTree(const size_type& __n, const std::vector<bool>& __rep) : _size(__n), _rep(__rep) { }
	BinaryWordPrefixFreeTree(const size_type& __n, const bool* __ptr) : _size(__n), _rep(__n,__ptr) { }
	BinaryWordPrefixFreeTree(const BinaryWordList& __list);
	
	/*!\brief The number of words in the list.
	 */
	size_type size() const { return _size; }

	/*!\brief Returns a constant forward iterator pointing to the beginning of the BinaryWordPrefixFreeTree.
	 */
	const_iterator begin() const;
	/*!\brief Returns a constant forward iterator pointing to the end of the BinaryWordPrefixFreeTree.
	 */
	const_iterator end() const;
      private:
	/* \brief The number of bytes needed to store the list.
	 */
	size_type storage_in_bytes() const { return _rep.size() / _bits_per_byte; }
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

    inline BinaryWordPrefixFreeTree::const_iterator BinaryWordPrefixFreeTree::begin() const { return ++(const_iterator(_rep.begin())); }
    inline BinaryWordPrefixFreeTree::const_iterator BinaryWordPrefixFreeTree::end() const { return const_iterator(_rep.end()); }



    BinaryWordPrefixFreeTree 
    inline
    join(const BinaryWordPrefixFreeTree& __tree1, const BinaryWordPrefixFreeTree& __tree2) {
	BinaryWordPrefixFreeTree result;
	std::vector<bool>& tree=result._rep;
	BinaryWord word;
	/*
	BinaryWordPrefixFreeTree::const_iterator iter1=__tree1.begin();
	BinaryWordPrefixFreeTree::const_iterator iter2=__tree2.begin();
	*/	
	std::vector<bool>::const_iterator iter1=__tree1._rep.begin();
	std::vector<bool>::const_iterator iter2=__tree2._rep.begin();

	static const bool _leaf = BinaryWordPrefixFreeTree::_leaf;
	static const bool _branch = BinaryWordPrefixFreeTree::_branch;
	static const bool _empty = BinaryWordPrefixFreeTree::_empty;
	static const bool _full =  BinaryWordPrefixFreeTree::_full;
	
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
    }
		
    /* Alternative representation of a binary word using STL bitset. */
    class BinaryWord_Bitset {
	friend std::ostream& operator<<(std::ostream&, const BinaryWord_Bitset&);
	typedef size_t size_type;
      private:
	size_type _size;
	std::bitset<256> _bits;
	static const size_type _max_size=255;
      public:
	BinaryWord_Bitset() : _size(0) { };
	BinaryWord_Bitset(size_type n, const void* p) : _size(n), _bits(*static_cast<const std::bitset<256>*>(p)) { }
	BinaryWord_Bitset(size_type n, const std::bitset<256>& bs) : _size(n), _bits(bs) { }
	BinaryWord_Bitset(const BinaryWord_Bitset& b) : _size(b._size), _bits(b._bits) { }
	
	size_type size() const { return _size; }
	bool operator[] (size_type n) const { return _bits[n]; }
	bool at(size_type n) const { assert(n<_size); return operator[](n); }
	const void* pointer() const { return &_bits; }    
	
	bool is_subword(const BinaryWord_Bitset& b) const {
	    if(size()>b.size()) { return false; }
	    for(size_type i=0; i!=size(); ++i) { if((*this)[i] != b[i]) { return false; } }
	    return true;
	}
    };
    
    inline
    std::ostream& operator<<(std::ostream& os, const BinaryWord& b) 
    {    
	os << b.size() << " :"; for(int i=0; i!=b.size(); ++i) { if(i%8==0) { os << " "; } os << b[i];  }
	return os;
    }
    
    inline
    std::ostream& operator<<(std::ostream& os, const BinaryWord_Bitset& b) 
    {    
	os << b.size() << " :"; for(int i=0; i!=b.size(); ++i) { if(i%8==0) { os << " "; } os << b[i];  }
	return os;
    }
    
}

#endif /* _BINARY_WORD_H */
