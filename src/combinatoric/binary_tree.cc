/***************************************************************************
 *            binary_tree.cc
 *
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

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

#include "combinatoric/binary_tree.h"

#include "base/array.h"
#include "base/stlio.h"

namespace Ariadne {
  namespace Combinatoric {

  inline size_type pow2(size_type n) { return 1<<n; }
  
  BinaryTree::BinaryTree(const Base::array<bool>& t)
    : _array(t.begin(),t.end()) 
  { 
    check(); 
  }
  
  BinaryTree::BinaryTree(const std::vector<bool>& t)
    : _array(t.begin(),t.end()) 
  { 
    check(); 
  }
  
  BinaryTree::BinaryTree(const size_type& n, const bool* ptr)
    : _array(n,ptr) 
  { 
      check();
  }
  
  BinaryTree::BinaryTree()
     : _array(1) { _array[0]=leaf; }

    BinaryTree::BinaryTree(size_type d)
      : _array(pow2(d+1)-1)
    {
      size_type n=_array.size();
      _array[n-1]=leaf;

      for(size_type i=0; i!=d; ++i) {
        size_type m=pow2(i+1)-1;
        for(size_type j=n-m; j!=n; ++j) {
          _array[j-m]=_array[j];
        }
        _array[n-2*m-1]=branch;
      }
    }

    void 
    BinaryTree::check() const { 
      size_type ends=1;
      for(BooleanArray::const_iterator iter=_array.begin(); iter!=_array.end(); ++iter) {
        /* Test for too many leaves */
        assert(ends!=0); 
        if(*iter==branch) {
          ends+=1;
        }
        else {
          ends-=1;
        }
      }
      assert(ends==0);
    }
    
    size_type 
    BinaryTree::depth() const { 
      size_type depth=0;
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        if(iter->size() > depth) {
          depth=iter->size();
        }
      }
      return depth;
    }
    
    
    void BinaryTree::advance(BinaryWord& word)
    {
      while(!word.empty() && word.back()==BinaryTree::right) {
        word.pop_back();
      }
      if(!word.empty()) {
        word.set_back(BinaryTree::right);
      }
    }
    
    
    
    BinaryTreeIterator::BinaryTreeIterator(BooleanArray::const_iterator i, BooleanArray::const_iterator e)
      : _position(i), _word()
    {
      if(i!=e) {
        while(*_position==BinaryTree::branch) { 
          _word.push_back(BinaryTree::left); 
          ++_position; 
        }
      }
    }
    
    void
    BinaryTreeIterator::increment() { 
      while( !_word.empty() && (_word.back() == BinaryTree::right) ) {
        _word.pop_back();
      }
      
      if( _word.empty() ) {
        ++_position;
        return;
      }
      
      _word.set_back(BinaryTree::right);
      ++_position;
      
      while( (*_position) == BinaryTree::branch ) {
        _word.push_back(BinaryTree::left);
        ++_position;
      }
      return;
    }
    
    
    void 
    BinaryTreeIterator::skip_subtree() 
    {
      size_type __len=_word.size();
      do {
        if(*_position == BinaryTree::branch) {
          _word.push_back(BinaryTree::left);
          ++_position;
        }
        else {
          ++_position;
          ++_position;
          if(_word.size()==__len) {
            return;
          }
          while(_word.back()==BinaryTree::right) {
            _word.pop_back();
            if(_word.size()==__len) {
              return;
            }
          }
          _word.set_back(BinaryTree::right);
        }
      }
      while(true);
    }
    
   
    std::ostream& 
    operator<<(std::ostream& os, const BinaryTree& t) {
      os << "BinaryTree( words=";
      Base::write_sequence(os, t.begin(), t.end());
      return os << ", sequence=" << t._array << " )"; 
    }
    
  }
}
