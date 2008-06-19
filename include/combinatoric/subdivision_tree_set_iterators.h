/***************************************************************************
 *            subdivision_tree_set_iterators.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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

/*! \file subdivision_tree_set.h
 *  \brief Cuboidal partition trees on a unit cuboid.
 */

#ifndef ARIADNE_SUBDIVISION_TREE_SET_ITERATORS_H
#define ARIADNE_SUBDIVISION_TREE_SET_ITERATORS_H

#include "subdivision_tree_set.h"

#include "base/iterator.h"

namespace Ariadne {

    
    inline uint log2(uint n) {
      uint r=0;  while(n>1) { ++r; n>>=1; } return r; 
    }

 
    class SubdivisionCellListSetIterator 
      : public boost::iterator_facade<SubdivisionCellListSetIterator,
                                      SubdivisionCell,
                                      boost::random_access_traversal_tag,
                                      SubdivisionCell>
    {
      typedef std::vector<word_type>::const_iterator base_iter;
      friend class boost::iterator_core_access;
     public:
      SubdivisionCellListSetIterator(dimension_type dim, base_iter curr) 
        : _dimension(dim), _curr(curr) { }
     private:
      bool equal(const SubdivisionCellListSetIterator& other) const { 
        return this->_curr==other._curr; }
      void increment() { ++_curr; }
      void decrement() { --_curr; } 
      void advance(difference_type n) { _curr+=n; } 
      difference_type distance_to(const SubdivisionCellListSetIterator& other) {
        return other._curr - this->_curr; }
      SubdivisionCell dereference() const { 
        return SubdivisionCell(this->_dimension,*this->_curr); }
     private:
      uint _dimension;
      base_iter _curr;
    };



    class SubdivisionMaskSetIterator 
      : public boost::iterator_facade<SubdivisionMaskSetIterator,
                                      SubdivisionCell,
                                      boost::forward_traversal_tag,
                                      SubdivisionCell>
    {
      typedef array<bool>::const_iterator base_iter;
      friend class boost::iterator_core_access;
     public:
      SubdivisionMaskSetIterator(dimension_type dim, base_iter curr, base_iter begin, base_iter end) 
        : _dimension(dim), _depth(log2(end-begin)) , _curr(curr), _begin(begin), _end(end)
      { while(_curr!=_end && !*_curr) { ++_curr; } }
     private:
      bool equal(const SubdivisionMaskSetIterator& other) const { 
        return this->_curr==other._curr; }
      void increment() { 
        do { ++_curr; } 
        while(_curr<_end && !*_curr); }
      void decrement() { 
        do { --_curr; } 
        while(_curr>=_begin && !*_curr); }
      SubdivisionCell dereference() const { 
        return SubdivisionCell(this->_dimension,this->_depth,this->_curr-this->_begin); }
     private:
      uint _dimension;
      uint _depth;
      base_iter _curr;
      base_iter _begin;
      base_iter _end;
     };





    class SubdivisionTreeSetIterator 
      : public boost::iterator_adaptor<SubdivisionTreeSetIterator,
                                       mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator>,
                                       SubdivisionCell,
                                       boost::use_default,
                                       SubdivisionCell>
      
    {
      typedef mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator> Base;
     public:
      SubdivisionTreeSetIterator(const SubdivisionSequence& ss,
                                   BinaryTree::const_iterator ti, 
                                   BooleanArray::const_iterator mi, 
                                   BooleanArray::const_iterator me) 
      //        : SubdivisionTreeSetIterator::iterator_adaptor_(mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator>(ti,mi,me)),
        : SubdivisionTreeSetIterator::iterator_adaptor_(Base(ti,mi,me)),
          _subdivisions(ss)
      { }
     private:
      friend class boost::iterator_core_access;
      SubdivisionCell dereference() const { 
        return SubdivisionCell(_subdivisions,*this->base_reference()); 
      }
     private:
      SubdivisionSequence _subdivisions;
    };
  
  
  
} // namespace Ariadne


#endif /* ARIADNE_SUBDIVISION_TREE_SET_ITERATORS_H */
