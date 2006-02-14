/***************************************************************************
 *            binary_tree.cc
 *
 *  16 November 2005
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

#include <base/binary_tree.h>

namespace Ariadne {
  BinaryTreeIterator& 
  BinaryTreeIterator::operator++() { 
    while( !_word.empty() && (_word.back() == right) ) {
      _word.pop_back();
    }
    
    if( _word.empty() ) {
      ++_position;
      return *this;
    }

    _word.set_back(right);
    ++_position;

    while( (*_position) == branch ) {
      _word.push_back(left);
      ++_position;
    }
    return *this;
  }
    
     
  void 
  BinaryTreeIterator::skip_subtree() 
  {
    size_type __len=_word.size();
    do {
      if(*_position == branch) {
        _word.push_back(0);
        ++_position;
      }
      else {
        ++_position;
        ++_position;
        if(_word.size()==__len) {
          return;
        }
        while(_word.back()==1) {
          _word.pop_back();
          if(_word.size()==__len) {
            return;
          }
        }
        _word.set_back(1);
      }
    }
    while(true);
  }

   
  std::ostream& 
  operator<<(std::ostream& os, const BinaryTree& t) {
    os << "BinaryTree( words=";
    write_sequence(os, t.begin(), t.end());
    return os << ", sequence=" << t._array << " )"; 
  }

}
